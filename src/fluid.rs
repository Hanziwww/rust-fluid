use rayon::prelude::*;

const U_FIELD: usize = 0;
const V_FIELD: usize = 1;
const S_FIELD: usize = 2;

// Standalone field sampling function (bilinear interpolation - optimized)
#[inline(always)]
fn sample_field(
    x: f32,
    y: f32,
    field: usize,
    field_data: &[f32],
    u_data: &[f32],
    v_data: &[f32],
    m_data: &[f32],
    num_x: usize,
    num_y: usize,
    h: f32,
) -> f32 {
    let n = num_y;
    let h1 = 1.0 / h;
    let h2 = 0.5 * h;

    let x = x.clamp(h, num_x as f32 * h);
    let y = y.clamp(h, num_y as f32 * h);

    let (dx, dy, f) = match field {
        U_FIELD => (0.0, h2, u_data),
        V_FIELD => (h2, 0.0, v_data),
        S_FIELD => (h2, h2, m_data),
        _ => (0.0, 0.0, field_data),
    };

    let x_adj = (x - dx) * h1;
    let y_adj = (y - dy) * h1;
    
    let x0 = x_adj.floor().min((num_x - 1) as f32) as usize;
    let y0 = y_adj.floor().min((num_y - 1) as f32) as usize;
    
    let tx = x_adj - x0 as f32;
    let ty = y_adj - y0 as f32;
    
    let x1 = (x0 + 1).min(num_x - 1);
    let y1 = (y0 + 1).min(num_y - 1);

    // Compute weights
    let sx = 1.0 - tx;
    let sy = 1.0 - ty;
    
    // Bilinear interpolation using FMA-friendly operations
    let v00 = f[x0 * n + y0];
    let v10 = f[x1 * n + y0];
    let v01 = f[x0 * n + y1];
    let v11 = f[x1 * n + y1];
    
    // Optimize: combine multiplications
    sx * sy * v00 + tx * sy * v10 + tx * ty * v11 + sx * ty * v01
}

#[derive(Clone)]
pub struct Fluid {
    pub density: f32,
    pub num_x: usize,
    pub num_y: usize,
    pub num_cells: usize,
    pub h: f32,
    pub u: Vec<f32>,
    pub v: Vec<f32>,
    pub p: Vec<f32>,
    pub s: Vec<f32>,
    pub m: Vec<f32>,
}

impl Fluid {
    pub fn new(density: f32, num_x: usize, num_y: usize, h: f32) -> Self {
        let num_x = num_x + 2;
        let num_y = num_y + 2;
        let num_cells = num_x * num_y;

        Fluid {
            density,
            num_x,
            num_y,
            num_cells,
            h,
            u: vec![0.0; num_cells],
            v: vec![0.0; num_cells],
            p: vec![0.0; num_cells],
            s: vec![0.0; num_cells],
            m: vec![1.0; num_cells],
        }
    }

    // Gravity integration (optimized - direct parallel modification)
    pub fn integrate(&mut self, dt: f32, gravity: f32) {
        let n = self.num_y;
        let gravity_dt = gravity * dt;

        // Direct parallel modification using indexed access
        self.v.par_iter_mut().enumerate().for_each(|(idx, v_val)| {
            let i = idx / n;
            let j = idx % n;
            
            if i >= 1 && i < self.num_x && j >= 1 && j < self.num_y - 1 {
                if self.s[idx] != 0.0 && self.s[i * n + j - 1] != 0.0 {
                    *v_val += gravity_dt;
                }
            }
        });
    }

    // Solve incompressibility (optimized with less allocations)
    pub fn solve_incompressibility(&mut self, num_iters: usize, dt: f32, over_relaxation: f32) {
        let n = self.num_y;
        let cp = self.density * self.h / dt;

        for _iter in 0..num_iters {
            // Use red-black checkerboard pattern for parallelization to avoid race conditions
            // Red cells (i+j is even)
            self.solve_incompressibility_pass(n, cp, over_relaxation, true);
            // Black cells (i+j is odd)
            self.solve_incompressibility_pass(n, cp, over_relaxation, false);
        }
    }

    fn solve_incompressibility_pass(&mut self, n: usize, cp: f32, over_relaxation: f32, red: bool) {
        use std::sync::Mutex;
        
        // Pre-allocate a buffer for updates to avoid repeated allocations
        let updates = Mutex::new(Vec::with_capacity(self.num_cells / 4));
        
        // Parallel computation of pressure updates
        (1..self.num_x - 1).into_par_iter().for_each(|i| {
            let mut local_updates = Vec::with_capacity(self.num_y);
            
            for j in 1..self.num_y - 1 {
                // Red-black checkerboard pattern
                if ((i + j) % 2 == 0) != red {
                    continue;
                }

                if self.s[i * n + j] == 0.0 {
                    continue;
                }

                let sx0 = self.s[(i - 1) * n + j];
                let sx1 = self.s[(i + 1) * n + j];
                let sy0 = self.s[i * n + j - 1];
                let sy1 = self.s[i * n + j + 1];
                let s_sum = sx0 + sx1 + sy0 + sy1;
                
                if s_sum == 0.0 {
                    continue;
                }

                let div = self.u[(i + 1) * n + j] - self.u[i * n + j]
                    + self.v[i * n + j + 1] - self.v[i * n + j];

                let p = -div / s_sum * over_relaxation;

                local_updates.push((i, j, p, sx0, sx1, sy0, sy1));
            }
            
            if !local_updates.is_empty() {
                updates.lock().unwrap().extend(local_updates);
            }
        });

        // Apply updates (serial - guaranteed correct order)
        let updates = updates.into_inner().unwrap();
        for (i, j, p, sx0, sx1, sy0, sy1) in updates {
            self.p[i * n + j] += cp * p;
            self.u[i * n + j] -= sx0 * p;
            self.u[(i + 1) * n + j] += sx1 * p;
            self.v[i * n + j] -= sy0 * p;
            self.v[i * n + j + 1] += sy1 * p;
        }
    }

    // Extrapolate boundary conditions
    pub fn extrapolate(&mut self) {
        let n = self.num_y;
        
        for i in 0..self.num_x {
            self.u[i * n] = self.u[i * n + 1];
            self.u[i * n + self.num_y - 1] = self.u[i * n + self.num_y - 2];
        }
        
        for j in 0..self.num_y {
            self.v[j] = self.v[n + j];
            self.v[(self.num_x - 1) * n + j] = self.v[(self.num_x - 2) * n + j];
        }
    }



    // Velocity advection (optimized - reduced allocations)
    pub fn advect_vel(&mut self, dt: f32) {
        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;
        let num_x = self.num_x;
        let num_y = self.num_y;

        // Clone data for parallel read-only access
        let u_clone = self.u.clone();
        let v_clone = self.v.clone();
        let s_clone = self.s.clone();
        let m_clone = self.m.clone();
        
        let new_u_vals: Vec<f32> = (0..self.num_cells).into_par_iter().map(|idx| {
            let i = idx / n;
            let j = idx % n;
            
            if i >= 1 && i < num_x && j >= 1 && j < num_y - 1 {
                if s_clone[i * n + j] != 0.0 && s_clone[(i - 1) * n + j] != 0.0 {
                    let x = i as f32 * h;
                    let y = j as f32 * h + h2;
                    let u = u_clone[i * n + j];
                    // avg_v inline
                    let v = (v_clone[(i - 1) * n + j] + v_clone[i * n + j]
                        + v_clone[(i - 1) * n + j + 1] + v_clone[i * n + j + 1]) * 0.25;
                    let new_x = x - dt * u;
                    let new_y = y - dt * v;
                    sample_field(new_x, new_y, U_FIELD, &u_clone, &u_clone, &v_clone, &m_clone, num_x, num_y, h)
                } else {
                    u_clone[idx]
                }
            } else {
                u_clone[idx]
            }
        }).collect();
        
        let new_v_vals: Vec<f32> = (0..self.num_cells).into_par_iter().map(|idx| {
            let i = idx / n;
            let j = idx % n;
            
            if i >= 1 && i < num_x - 1 && j >= 1 && j < num_y {
                if s_clone[i * n + j] != 0.0 && s_clone[i * n + j - 1] != 0.0 {
                    let x = i as f32 * h + h2;
                    let y = j as f32 * h;
                    // avg_u inline
                    let u = (u_clone[i * n + j - 1] + u_clone[i * n + j]
                        + u_clone[(i + 1) * n + j - 1] + u_clone[(i + 1) * n + j]) * 0.25;
                    let v = v_clone[i * n + j];
                    let new_x = x - dt * u;
                    let new_y = y - dt * v;
                    sample_field(new_x, new_y, V_FIELD, &v_clone, &u_clone, &v_clone, &m_clone, num_x, num_y, h)
                } else {
                    v_clone[idx]
                }
            } else {
                v_clone[idx]
            }
        }).collect();

        self.u = new_u_vals;
        self.v = new_v_vals;
    }

    // Smoke advection (optimized - reduced allocations)
    pub fn advect_smoke(&mut self, dt: f32) {
        let n = self.num_y;
        let h = self.h;
        let h2 = 0.5 * h;
        let num_x = self.num_x;
        let num_y = self.num_y;

        // Clone data for parallel read-only access
        let u_clone = self.u.clone();
        let v_clone = self.v.clone();
        let s_clone = self.s.clone();
        let m_clone = self.m.clone();

        let new_m_vals: Vec<f32> = (0..self.num_cells).into_par_iter().map(|idx| {
            let i = idx / n;
            let j = idx % n;
            
            if i >= 1 && i < num_x - 1 && j >= 1 && j < num_y - 1 {
                if s_clone[i * n + j] != 0.0 {
                    let u = (u_clone[i * n + j] + u_clone[(i + 1) * n + j]) * 0.5;
                    let v = (v_clone[i * n + j] + v_clone[i * n + j + 1]) * 0.5;
                    let x = i as f32 * h + h2 - dt * u;
                    let y = j as f32 * h + h2 - dt * v;
                    sample_field(x, y, S_FIELD, &m_clone, &u_clone, &v_clone, &m_clone, num_x, num_y, h)
                } else {
                    m_clone[idx]
                }
            } else {
                m_clone[idx]
            }
        }).collect();

        self.m = new_m_vals;
    }

    // Main simulation step
    pub fn simulate(&mut self, dt: f32, gravity: f32, num_iters: usize, over_relaxation: f32) {
        self.integrate(dt, gravity);
        self.p.fill(0.0);
        self.solve_incompressibility(num_iters, dt, over_relaxation);
        self.extrapolate();
        self.advect_vel(dt);
        self.advect_smoke(dt);
    }

    // Set obstacle (optimized - parallel processing)
    pub fn set_obstacle(&mut self, x: f32, y: f32, radius: f32, _reset: bool, vx: f32, vy: f32, scene_nr: usize, frame_nr: usize) {
        let n = self.num_y;
        let r = radius;
        let h = self.h;
        let num_x = self.num_x;
        let num_y = self.num_y;
        let r_squared = r * r;
        let sin_val = 0.5 + 0.5 * ((0.1 * frame_nr as f32).sin());

        // Collect updates in parallel
        use std::sync::Mutex;
        let updates = Mutex::new(Vec::new());

        (1..num_x - 2).into_par_iter().for_each(|i| {
            let mut local_updates = Vec::new();
            
            for j in 1..num_y - 2 {
                let idx = i * n + j;
                let dx = (i as f32 + 0.5) * h - x;
                let dy = (j as f32 + 0.5) * h - y;

                if dx * dx + dy * dy < r_squared {
                    let m_val = if scene_nr == 2 { sin_val } else { 1.0 };
                    local_updates.push((idx, m_val));
                }
            }
            
            if !local_updates.is_empty() {
                updates.lock().unwrap().extend(local_updates);
            }
        });

        // First, reset all s values in the working area
        for i in 1..num_x - 2 {
            for j in 1..num_y - 2 {
                self.s[i * n + j] = 1.0;
            }
        }

        // Apply updates serially (necessary for correctness)
        let updates = updates.into_inner().unwrap();
        for (idx, m_val) in updates {
            let i = idx / n;
            let j = idx % n;
            
            self.s[idx] = 0.0;
            self.m[idx] = m_val;
            self.u[idx] = vx;
            self.u[(i + 1) * n + j] = vx;
            self.v[idx] = vy;
            self.v[i * n + j + 1] = vy;
        }
    }
}

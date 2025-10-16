use crate::fluid::Fluid;
use crate::profiler::Profiler;

#[derive(Clone, Copy, PartialEq)]
pub enum SceneType {
    Tank = 0,
    WindTunnel = 1,
    Paint = 2,
    HiresTunnel = 3,
}

pub struct Scene {
    pub gravity: f32,
    pub dt: f32,
    pub num_iters: usize,
    pub frame_nr: usize,
    pub over_relaxation: f32,
    pub obstacle_x: f32,
    pub obstacle_y: f32,
    pub obstacle_radius: f32,
    pub paused: bool,
    pub scene_nr: SceneType,
    pub show_obstacle: bool,
    pub show_streamlines: bool,
    pub show_velocities: bool,
    pub show_pressure: bool,
    pub show_smoke: bool,
    pub fluid: Option<Fluid>,
    pub sim_height: f32,
    pub sim_width: f32,
    pub profiler: Option<Profiler>,
}

impl Default for Scene {
    fn default() -> Self {
        Scene {
            gravity: -9.81,
            dt: 1.0 / 120.0,
            num_iters: 100,
            frame_nr: 0,
            over_relaxation: 1.9,
            obstacle_x: 0.0,
            obstacle_y: 0.0,
            obstacle_radius: 0.15,
            paused: false,
            scene_nr: SceneType::WindTunnel,
            show_obstacle: false,
            show_streamlines: false,
            show_velocities: false,
            show_pressure: false,
            show_smoke: true,
            fluid: None,
            sim_height: 1.1,
            sim_width: 1.0,
            profiler: Some(Profiler::new()),
        }
    }
}

impl Scene {
    pub fn new(width: f32, height: f32) -> Self {
        let sim_height = 1.1;
        let sim_width = width / height * sim_height;
        
        Scene {
            sim_height,
            sim_width,
            ..Default::default()
        }
    }

    pub fn setup(&mut self, scene_nr: SceneType) {
        self.scene_nr = scene_nr;
        self.obstacle_radius = 0.15;
        self.over_relaxation = 1.9;
        self.dt = 1.0 / 60.0;
        self.num_iters = 40;

        let res = match scene_nr {
            SceneType::Tank => 200,
            SceneType::HiresTunnel => 200,
            _ => 100,
        };

        let domain_height = 1.0;
        let domain_width = domain_height / self.sim_height * self.sim_width;
        let h = domain_height / res as f32;

        let num_x = (domain_width / h).floor() as usize;
        let num_y = (domain_height / h).floor() as usize;

        let density = 1000.0;

        let mut f = Fluid::new(density, num_x, num_y, h);
        let n = f.num_y;

        match scene_nr {
            SceneType::Tank => {
                // Tank scene
                for i in 0..f.num_x {
                    for j in 0..f.num_y {
                        let s = if i == 0 || i == f.num_x - 1 || j == 0 {
                            0.0 // Solid
                        } else {
                            1.0 // Fluid
                        };
                        f.s[i * n + j] = s;
                    }
                }
                self.gravity = -9.81;
                self.show_pressure = true;
                self.show_smoke = false;
                self.show_streamlines = false;
                self.show_velocities = false;
            }

            SceneType::WindTunnel | SceneType::HiresTunnel => {
                // Wind tunnel scene
                let in_vel = 2.0;
                for i in 0..f.num_x {
                    for j in 0..f.num_y {
                        let s = if i == 0 || j == 0 || j == f.num_y - 1 {
                            0.0 // Solid
                        } else {
                            1.0 // Fluid
                        };
                        f.s[i * n + j] = s;

                        if i == 1 {
                            f.u[i * n + j] = in_vel;
                        }
                    }
                }

                let pipe_h = (0.1 * f.num_y as f32) as usize;
                let min_j = (0.5 * f.num_y as f32 - 0.5 * pipe_h as f32).floor() as usize;
                let max_j = (0.5 * f.num_y as f32 + 0.5 * pipe_h as f32).floor() as usize;

                for j in min_j..max_j {
                    f.m[j] = 0.0;
                }

                self.set_obstacle(0.4, 0.5, true, 0.0, 0.0);

                self.gravity = 0.0;
                self.show_pressure = false;
                self.show_smoke = true;
                self.show_streamlines = false;
                self.show_velocities = false;

                if scene_nr == SceneType::HiresTunnel {
                    self.dt = 1.0 / 120.0;
                    self.num_iters = 100;
                    self.show_pressure = true;
                }
            }

            SceneType::Paint => {
                // Paint scene
                self.gravity = 0.0;
                self.over_relaxation = 1.0;
                self.show_pressure = false;
                self.show_smoke = true;
                self.show_streamlines = false;
                self.show_velocities = false;
                self.obstacle_radius = 0.1;
            }
        }

        self.fluid = Some(f);
    }

    pub fn set_obstacle(&mut self, x: f32, y: f32, reset: bool, vx_in: f32, vy_in: f32) {
        if let Some(ref mut fluid) = self.fluid {
            let (vx, vy) = if reset {
                (0.0, 0.0)
            } else {
                let vx = (x - self.obstacle_x) / self.dt;
                let vy = (y - self.obstacle_y) / self.dt;
                (vx, vy)
            };

            self.obstacle_x = x;
            self.obstacle_y = y;

            fluid.set_obstacle(
                x,
                y,
                self.obstacle_radius,
                reset,
                vx + vx_in,
                vy + vy_in,
                self.scene_nr as usize,
                self.frame_nr,
            );

            self.show_obstacle = true;
        }
    }

    pub fn simulate(&mut self) {
        if !self.paused {
            if let Some(ref mut fluid) = self.fluid {
                if let Some(ref mut profiler) = self.profiler {
                    profiler.start("integrate");
                    fluid.integrate(self.dt, self.gravity);
                    profiler.end();
                    
                    fluid.p.fill(0.0);
                    
                    profiler.start("incompressibility");
                    fluid.solve_incompressibility(self.num_iters, self.dt, self.over_relaxation);
                    profiler.end();
                    
                    profiler.start("extrapolate");
                    fluid.extrapolate();
                    profiler.end();
                    
                    profiler.start("advect_vel");
                    fluid.advect_vel(self.dt);
                    profiler.end();
                    
                    profiler.start("advect_smoke");
                    fluid.advect_smoke(self.dt);
                    profiler.end();
                } else {
                    fluid.simulate(self.dt, self.gravity, self.num_iters, self.over_relaxation);
                }
                self.frame_nr += 1;
            }
        }
    }
}

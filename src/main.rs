mod fluid;
mod scene;
mod profiler;

use eframe::egui;
use scene::{Scene, SceneType};
use profiler::Profiler;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 800.0])
            .with_title("Fluid Simulator - Rust Parallel Edition"),
        ..Default::default()
    };

    eframe::run_native(
        "Fluid Simulator",
        options,
        Box::new(|cc| Ok(Box::new(FluidSimApp::new(cc)))),
    )
}

struct FluidSimApp {
    scene: Scene,
    mouse_pos: Option<egui::Pos2>,
    dragging: bool,
    prev_mouse_pos: Option<egui::Pos2>,
    profiler: Profiler,
    show_profiler: bool,
}

impl FluidSimApp {
    fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        let mut app = FluidSimApp {
            scene: Scene::new(1200.0, 700.0),
            mouse_pos: None,
            dragging: false,
            prev_mouse_pos: None,
            profiler: Profiler::new(),
            show_profiler: false,
        };
        app.scene.setup(SceneType::WindTunnel);
        app
    }

    fn draw_fluid(&self, painter: &egui::Painter, rect: egui::Rect) {
        use rayon::prelude::*;
        
        if let Some(ref fluid) = self.scene.fluid {
            let c_scale = rect.height() / self.scene.sim_height;
            let n = fluid.num_y;
            let h = fluid.h;
            let cell_scale = 1.1;

            // Calculate pressure range using parallel iterator
            let (min_p, max_p) = if self.scene.show_pressure {
                use rayon::prelude::*;
                let (min_p, max_p) = fluid.p.par_iter()
                    .fold(|| (f32::INFINITY, f32::NEG_INFINITY), |acc, &p| {
                        (acc.0.min(p), acc.1.max(p))
                    })
                    .reduce(|| (f32::INFINITY, f32::NEG_INFINITY), |a, b| {
                        (a.0.min(b.0), a.1.max(b.1))
                    });
                (min_p, max_p)
            } else {
                (0.0, 1.0)
            };

            // Optimized: Single-level parallel iteration without nested par_iter
            let total_cells = fluid.num_x * fluid.num_y;
            let show_pressure = self.scene.show_pressure;
            let show_smoke = self.scene.show_smoke;
            let scene_nr = self.scene.scene_nr;
            
            let cells: Vec<_> = (0..total_cells).into_par_iter()
                .map(|idx| {
                    let i = idx / fluid.num_y;
                    let j = idx % fluid.num_y;
                    
                    let color = if show_pressure {
                        let p = fluid.p[i * n + j];
                        let s = fluid.m[i * n + j];
                        let mut rgb = get_sci_color(p, min_p, max_p);
                        
                        if show_smoke {
                            rgb[0] = (rgb[0] as f32 - 255.0 * s).max(0.0) as u8;
                            rgb[1] = (rgb[1] as f32 - 255.0 * s).max(0.0) as u8;
                            rgb[2] = (rgb[2] as f32 - 255.0 * s).max(0.0) as u8;
                        }
                        egui::Color32::from_rgb(rgb[0], rgb[1], rgb[2])
                    } else if show_smoke {
                        let s = fluid.m[i * n + j];
                        if scene_nr == SceneType::Paint {
                            let rgb = get_sci_color(s, 0.0, 1.0);
                            egui::Color32::from_rgb(rgb[0], rgb[1], rgb[2])
                        } else {
                            let val = (255.0 * s) as u8;
                            egui::Color32::from_rgb(val, val, val)
                        }
                    } else if fluid.s[i * n + j] == 0.0 {
                        egui::Color32::BLACK
                    } else {
                        egui::Color32::WHITE
                    };

                    let x = rect.left() + i as f32 * h * c_scale;
                    let y = rect.bottom() - (j as f32 + 1.0) * h * c_scale;
                    let cx = (c_scale * cell_scale * h) + 1.0;
                    let cy = (c_scale * cell_scale * h) + 1.0;

                    let cell_rect = egui::Rect::from_min_size(
                        egui::pos2(x, y),
                        egui::vec2(cx, cy),
                    );
                    (cell_rect, color)
                })
                .collect();

            // Batch draw all cells
            for (cell_rect, color) in cells {
                painter.rect_filled(cell_rect, 0.0, color);
            }

            // Draw velocity field
            if self.scene.show_velocities {
                let scale = 0.02;
                for i in 0..fluid.num_x {
                    for j in 0..fluid.num_y {
                        let u = fluid.u[i * n + j];
                        let v = fluid.v[i * n + j];

                        let x = rect.left() + i as f32 * h * c_scale;
                        let y = rect.bottom() - (j as f32 + 0.5) * h * c_scale;

                        // Draw u component
                        let x0 = x;
                        let x1 = x + u * scale * c_scale;
                        painter.line_segment(
                            [egui::pos2(x0, y), egui::pos2(x1, y)],
                            egui::Stroke::new(1.0, egui::Color32::BLACK),
                        );

                        // Draw v component
                        let x_mid = rect.left() + (i as f32 + 0.5) * h * c_scale;
                        let y0 = rect.bottom() - j as f32 * h * c_scale;
                        let y1 = y0 - v * scale * c_scale;
                        painter.line_segment(
                            [egui::pos2(x_mid, y0), egui::pos2(x_mid, y1)],
                            egui::Stroke::new(1.0, egui::Color32::BLACK),
                        );
                    }
                }
            }

            // Draw streamlines
            if self.scene.show_streamlines {
                let _seg_len = h * 0.2;
                let num_segs = 15;

                for i in (1..fluid.num_x - 1).step_by(5) {
                    for j in (1..fluid.num_y - 1).step_by(5) {
                        let mut x = (i as f32 + 0.5) * h;
                        let mut y = (j as f32 + 0.5) * h;

                        let mut points = Vec::new();
                        points.push(egui::pos2(
                            rect.left() + x * c_scale,
                            rect.bottom() - y * c_scale,
                        ));

                        for _ in 0..num_segs {
                            let u = fluid_sample_field(fluid, x, y, 0);
                            let v = fluid_sample_field(fluid, x, y, 1);
                            x += u * 0.01;
                            y += v * 0.01;
                            if x > fluid.num_x as f32 * h {
                                break;
                            }
                            points.push(egui::pos2(
                                rect.left() + x * c_scale,
                                rect.bottom() - y * c_scale,
                            ));
                        }

                        painter.add(egui::Shape::line(
                            points,
                            egui::Stroke::new(1.0, egui::Color32::BLACK),
                        ));
                    }
                }
            }

            // Draw obstacle
            if self.scene.show_obstacle {
                let r = (self.scene.obstacle_radius + h) * c_scale;
                let center = egui::pos2(
                    rect.left() + self.scene.obstacle_x * c_scale,
                    rect.bottom() - self.scene.obstacle_y * c_scale,
                );

                let color = if self.scene.show_pressure {
                    egui::Color32::BLACK
                } else {
                    egui::Color32::from_rgb(221, 221, 221)
                };

                painter.circle_filled(center, r, color);
                painter.circle_stroke(center, r, egui::Stroke::new(3.0, egui::Color32::BLACK));
            }

            // Display pressure information
            if self.scene.show_pressure {
                let text = format!("Pressure: {:.0} - {:.0} N/m", min_p, max_p);
                painter.text(
                    egui::pos2(rect.left() + 10.0, rect.top() + 10.0),
                    egui::Align2::LEFT_TOP,
                    text,
                    egui::FontId::proportional(16.0),
                    egui::Color32::BLACK,
                );
            }
        }
    }
}

impl eframe::App for FluidSimApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("Control Panel").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if ui.button("Wind Tunnel").clicked() {
                    self.scene.setup(SceneType::WindTunnel);
                }
                if ui.button("Hires Tunnel").clicked() {
                    self.scene.setup(SceneType::HiresTunnel);
                }
                if ui.button("Tank").clicked() {
                    self.scene.setup(SceneType::Tank);
                }
                if ui.button("Paint").clicked() {
                    self.scene.setup(SceneType::Paint);
                }

                ui.separator();

                ui.checkbox(&mut self.scene.show_streamlines, "Streamlines");
                ui.checkbox(&mut self.scene.show_velocities, "Velocities");
                ui.checkbox(&mut self.scene.show_pressure, "Pressure");
                ui.checkbox(&mut self.scene.show_smoke, "Smoke");

                ui.separator();

                let over_relax = &mut self.scene.over_relaxation;
                let checked = *over_relax > 1.0;
                let mut new_checked = checked;
                ui.checkbox(&mut new_checked, "Overrelax");
                if new_checked != checked {
                    *over_relax = if new_checked { 1.9 } else { 1.0 };
                }

                ui.separator();

                if ui.button(if self.scene.paused { "▶ Resume" } else { "⏸ Pause" }).clicked() {
                    self.scene.paused = !self.scene.paused;
                }

                ui.separator();

                ui.checkbox(&mut self.show_profiler, "Show Profiler");

                ui.label(format!("Frame: {}", self.scene.frame_nr));
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let (response, painter) = ui.allocate_painter(
                ui.available_size(),
                egui::Sense::click_and_drag(),
            );

            let rect = response.rect;

            // Handle mouse interaction
            if let Some(pos) = response.interact_pointer_pos() {
                let c_scale = rect.height() / self.scene.sim_height;
                let x = (pos.x - rect.left()) / c_scale;
                let y = (rect.bottom() - pos.y) / c_scale;

                if response.drag_started() {
                    self.dragging = true;
                    self.prev_mouse_pos = Some(pos);
                    self.scene.set_obstacle(x, y, true, 0.0, 0.0);
                } else if self.dragging && response.dragged() {
                    if let Some(prev_pos) = self.prev_mouse_pos {
                        let _prev_x = (prev_pos.x - rect.left()) / c_scale;
                        let _prev_y = (rect.bottom() - prev_pos.y) / c_scale;
                        self.scene.set_obstacle(x, y, false, 0.0, 0.0);
                    }
                    self.prev_mouse_pos = Some(pos);
                }

                self.mouse_pos = Some(pos);
            }

            if response.drag_stopped() {
                self.dragging = false;
            }

            // Draw background
            painter.rect_filled(rect, 0.0, egui::Color32::WHITE);

            // Draw fluid
            self.profiler.start("draw_fluid");
            self.draw_fluid(&painter, rect);
            self.profiler.end();

            // Simulate
            self.profiler.start("simulate_total");
            self.scene.simulate();
            self.profiler.end();
        });

        // Show profiler if enabled
        if self.show_profiler {
            egui::Window::new("Performance Profiler")
                .default_width(600.0)
                .default_height(400.0)
                .show(ctx, |ui| {
                    ui.heading("Overall Performance");
                    ui.label(self.profiler.report());
                    
                    if let Some(ref scene_profiler) = self.scene.profiler {
                        ui.separator();
                        ui.heading("Simulation Details");
                        ui.label(scene_profiler.report());
                    }
                    
                    if ui.button("Reset Stats").clicked() {
                        self.profiler.reset();
                        if let Some(ref mut scene_profiler) = self.scene.profiler {
                            scene_profiler.reset();
                        }
                    }
                });
        }

        // Request continuous repaint
        ctx.request_repaint();
    }
}

// Helper function: Gnuplot2 colormap (matplotlib-style)
// Maps values from min to max using the classic gnuplot2 color scheme
// Color progression: dark blue -> purple -> red -> orange -> yellow -> light yellow
fn get_sci_color(val: f32, min_val: f32, max_val: f32) -> [u8; 3] {
    let val = val.max(min_val).min(max_val - 0.0001);
    let d = max_val - min_val;
    let t = if d == 0.0 { 0.5 } else { (val - min_val) / d };
    
    // Gnuplot2 colormap approximation (piecewise linear)
    let r = if t < 0.13 {
        0.0
    } else if t < 0.25 {
        (t - 0.13) / 0.12
    } else if t < 0.5 {
        1.0
    } else if t < 0.75 {
        1.0
    } else {
        1.0 - (t - 0.75) / 0.5
    };
    
    let g = if t < 0.13 {
        0.0
    } else if t < 0.38 {
        (t - 0.13) / 0.25 * 0.6
    } else if t < 0.5 {
        0.6 + (t - 0.38) / 0.12 * 0.4
    } else if t < 0.63 {
        1.0
    } else if t < 0.88 {
        1.0 - (t - 0.63) / 0.25 * 0.4
    } else {
        0.6 - (t - 0.88) / 0.12 * 0.6
    };
    
    let b = if t < 0.13 {
        0.54 + t / 0.13 * 0.46
    } else if t < 0.25 {
        1.0 - (t - 0.13) / 0.12
    } else if t < 0.5 {
        0.0
    } else if t < 0.75 {
        0.0
    } else {
        (t - 0.75) / 0.25
    };

    [
        (255.0 * r.clamp(0.0, 1.0)) as u8,
        (255.0 * g.clamp(0.0, 1.0)) as u8,
        (255.0 * b.clamp(0.0, 1.0)) as u8,
    ]
}

// Helper function: Field sampling
fn fluid_sample_field(fluid: &fluid::Fluid, x: f32, y: f32, field: usize) -> f32 {
    let n = fluid.num_y;
    let h = fluid.h;
    let h1 = 1.0 / h;
    let h2 = 0.5 * h;

    let x = x.max(h).min(fluid.num_x as f32 * h);
    let y = y.max(h).min(fluid.num_y as f32 * h);

    let (dx, dy, f) = match field {
        0 => (0.0, h2, &fluid.u),
        1 => (h2, 0.0, &fluid.v),
        2 => (h2, h2, &fluid.m),
        _ => (0.0, 0.0, &fluid.m),
    };

    let x0 = ((x - dx) * h1).floor().min((fluid.num_x - 1) as f32) as usize;
    let tx = ((x - dx) - x0 as f32 * h) * h1;
    let x1 = (x0 + 1).min(fluid.num_x - 1);

    let y0 = ((y - dy) * h1).floor().min((fluid.num_y - 1) as f32) as usize;
    let ty = ((y - dy) - y0 as f32 * h) * h1;
    let y1 = (y0 + 1).min(fluid.num_y - 1);

    let sx = 1.0 - tx;
    let sy = 1.0 - ty;

    sx * sy * f[x0 * n + y0]
        + tx * sy * f[x1 * n + y0]
        + tx * ty * f[x1 * n + y1]
        + sx * ty * f[x0 * n + y1]
}

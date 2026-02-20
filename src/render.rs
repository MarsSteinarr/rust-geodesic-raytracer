use anyhow::Context;
use image::RgbImage;
use rayon::prelude::*;

use crate::camera;
use crate::cli::Args;
use crate::cli::SkyMode;
use crate::geodesic;
use crate::integrator;
use crate::math;
use crate::metric;
use crate::sky;

fn build_camera(args: &Args, aspect: f64) -> camera::Camera {
    let m = args.m;

    let default_pos = [90.0 * m, 6.0 * m, 0.0 * m];
    let default_look = [0.0, 0.0, 0.0];
    let default_up = [0.0, 1.0, 0.0];

    let pos = args
        .cam_pos
        .as_ref()
        .and_then(|v| {
            if v.len() == 3 {
                Some([v[0] * m, v[1] * m, v[2] * m])
            } else {
                None
            }
        })
        .unwrap_or(default_pos);

    let look = args
        .cam_look
        .as_ref()
        .and_then(|v| {
            if v.len() == 3 {
                Some([v[0] * m, v[1] * m, v[2] * m])
            } else {
                None
            }
        })
        .unwrap_or(default_look);

    let up = args
        .cam_up
        .as_ref()
        .and_then(|v| {
            if v.len() == 3 {
                Some([v[0], v[1], v[2]])
            } else {
                None
            }
        })
        .unwrap_or(default_up);

    camera::make_camera(pos, look, up, args.fov_deg, aspect)
}

#[derive(Debug, Clone)]
pub struct DebugRayResult {
    pub steps: u32,
    pub event: String,
    pub max_abs_h: f64,
    pub final_r: f64,
}

#[derive(Clone, Copy, Debug)]
enum Event {
    Capture,
    Escape,
    Disk,
    MaxSteps,
}

fn shade_disk(r_xy: f64, r_in: f64, r_out: f64) -> [u8; 3] {
    let t = if r_out > r_in {
        ((r_xy - r_in) / (r_out - r_in)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let e = (1.0 - t).sqrt();

    let base = [255.0, 140.0, 40.0];
    let c0 = (base[0] * e).clamp(0.0, 255.0) as u8;
    let c1 = (base[1] * e).clamp(0.0, 255.0) as u8;
    let c2 = (base[2] * e).clamp(0.0, 255.0) as u8;
    [c0, c1, c2]
}

fn trace_ray_internal(
    px: u32,
    py: u32,
    cam: &camera::Camera,
    args: &Args,
) -> (Event, u32, f64, f64, [u8; 3]) {
    let m = args.m;
    let dir = camera::ray_dir(cam, px, py, args.width, args.height);
    let x0 = [0.0, cam.pos[0], cam.pos[1], cam.pos[2]];

    let k_up = [1.0, dir[0], dir[1], dir[2]];
    let p_flat = [-k_up[0], k_up[1], k_up[2], k_up[3]];

    let p0 = match geodesic::project_p0(x0, p_flat, m) {
        Some(p) => p,
        None => {
            return (Event::Capture, 0, 0.0, metric::r(x0), [0, 0, 0]);
        }
    };

    let mut y = [0.0; 8];
    y[0] = x0[0];
    y[1] = x0[1];
    y[2] = x0[2];
    y[3] = x0[3];
    y[4] = p0[0];
    y[5] = p0[1];
    y[6] = p0[2];
    y[7] = p0[3];

    let params = integrator::IntegratorParams {
        newton_iters: args.newton_iters,
        newton_tol: args.newton_tol,
    };

    let mut max_abs_h = 0.0;

    let disk_enabled = args.disk;

    for step in 0..args.max_steps {
        let x = [y[0], y[1], y[2], y[3]];
        let rr = metric::r(x);

        if rr < 1e-6 {
            return (Event::Capture, step, max_abs_h, rr, [0, 0, 0]);
        }
        if rr <= args.r_cap * m {
            return (Event::Capture, step, max_abs_h, rr, [0, 0, 0]);
        }
        if rr >= args.r_sky * m {
            let dy = geodesic::derivs(y, m);
            let dir3 = math::normalize3([dy[1], dy[2], dy[3]]);
            let col = match args.sky {
                SkyMode::Black => [0, 0, 0],
                SkyMode::Checker => sky::sample_checker(dir3, args.checker_squares),
                SkyMode::Stars => sky::sample_stars(dir3),
            };
            return (Event::Escape, step, max_abs_h, rr, col);
        }

        let y_next = match integrator::implicit_midpoint_step(y, args.step_h, m, params) {
            Some(v) => v,
            None => {
                return (Event::Capture, step, max_abs_h, rr, [0, 0, 0]);
            }
        };

        if disk_enabled {
            let h = args.disk_half_thickness * m;
            let x0p = y[1];
            let y0p = y[2];
            let z0p = y[3];
            let x1p = y_next[1];
            let y1p = y_next[2];
            let z1p = y_next[3];

            if h.is_finite()
                && z0p.is_finite()
                && z1p.is_finite()
                && x0p.is_finite()
                && y0p.is_finite()
                && x1p.is_finite()
                && y1p.is_finite()
            {
                let in0 = y0p.abs() <= h;
                let in1 = y1p.abs() <= h;
                let cross_pos = (y0p - h) * (y1p - h) <= 0.0;
                let cross_neg = (y0p + h) * (y1p + h) <= 0.0;
                let crosses_slab = in0 || in1 || cross_pos || cross_neg;

                if crosses_slab {
                    let dy = y1p - y0p;
                    if dy.abs() < 1e-12 && !(in0 || in1) {
                        continue;
                    }
                    let s = if (y0p * y1p <= 0.0) && dy.abs() > 0.0 {
                        y0p / (y0p - y1p)
                    } else if dy.abs() < 1e-12 {
                        if in0 { 0.0 } else { 1.0 }
                    } else if dy.abs() > 0.0 {
                        if cross_pos && !cross_neg {
                            (h - y0p) / dy
                        } else if cross_neg && !cross_pos {
                            (-h - y0p) / dy
                        } else if in0 && in1 {
                            let a0 = y0p.abs();
                            let a1 = y1p.abs();
                            if a0 <= a1 { 0.0 } else { 1.0 }
                        } else {
                            0.5
                        }
                    } else {
                        0.5
                    };

                    let s = s.clamp(0.0, 1.0);
                    let xh = x0p + s * (x1p - x0p);
                    let zh = z0p + s * (z1p - z0p);
                    let rho = (xh * xh + zh * zh).sqrt();

                    let r_in = args.disk_r_in * m;
                    let r_out = args.disk_r_out * m;
                    if rho.is_finite() && rho >= r_in && rho <= r_out {
                        let col = shade_disk(rho, r_in, r_out);
                        return (Event::Disk, step, max_abs_h, rr, col);
                    }
                }
            }
        }

        let x_new = [y_next[0], y_next[1], y_next[2], y_next[3]];
        let rr_new = metric::r(x_new);
        if rr_new < 1e-6 {
            return (Event::Capture, step, max_abs_h, rr_new, [0, 0, 0]);
        }
        if rr_new <= args.r_cap * m {
            return (Event::Capture, step, max_abs_h, rr_new, [0, 0, 0]);
        }
        let mut p_new = [y_next[4], y_next[5], y_next[6], y_next[7]];
        p_new = match geodesic::project_p0(x_new, p_new, m) {
            Some(p) => p,
            None => {
                let rr2 = metric::r(x_new);
                return (Event::Capture, step, max_abs_h, rr2, [0, 0, 0]);
            }
        };

        y = [
            x_new[0], x_new[1], x_new[2], x_new[3], p_new[0], p_new[1], p_new[2], p_new[3],
        ];

        let hval = geodesic::hamiltonian(x_new, p_new, m).abs();
        if hval.is_finite() && hval > max_abs_h {
            max_abs_h = hval;
        }

        if !y.iter().all(|v| v.is_finite()) {
            let rr2 = metric::r(x_new);
            return (Event::Capture, step, max_abs_h, rr2, [0, 0, 0]);
        }
    }

    let x = [y[0], y[1], y[2], y[3]];
    let rr = metric::r(x);
    (Event::MaxSteps, args.max_steps, max_abs_h, rr, [0, 0, 0])
}

pub fn trace_ray_debug(px: u32, py: u32, args: &Args) -> anyhow::Result<DebugRayResult> {
    if px >= args.width || py >= args.height {
        anyhow::bail!("debug ray pixel out of bounds");
    }

    let cam = build_camera(args, args.width as f64 / args.height as f64);
    let (ev, steps, max_abs_h, final_r, _) = trace_ray_internal(px, py, &cam, args);
    let event = match ev {
        Event::Capture => "capture",
        Event::Escape => "escape",
        Event::Disk => "disk",
        Event::MaxSteps => "max_steps",
    };

    Ok(DebugRayResult {
        steps,
        event: event.to_string(),
        max_abs_h,
        final_r,
    })
}

pub fn trace_ray(px: u32, py: u32, args: &Args) -> [u8; 3] {
    let cam = build_camera(args, args.width as f64 / args.height as f64);
    let (ev, _steps, _max_abs_h, _final_r, col) = trace_ray_internal(px, py, &cam, args);
    match ev {
        Event::Escape => col,
        Event::Disk => col,
        Event::Capture => [0, 0, 0],
        Event::MaxSteps => [0, 0, 0],
    }
}

pub fn render_image(args: &Args) -> anyhow::Result<RgbImage> {
    let width = args.width;
    let height = args.height;

    let cam = build_camera(args, width as f64 / height as f64);

    let mut buf = vec![0u8; (width as usize) * (height as usize) * 3];

    buf.par_chunks_mut(3)
        .enumerate()
        .for_each(|(idx, pix)| {
            let x = (idx as u32) % width;
            let y = (idx as u32) / width;
            let (ev, _steps, _max_abs_h, _final_r, col) = trace_ray_internal(x, y, &cam, args);
            let c = match ev {
                Event::Escape => col,
                Event::Disk => col,
                Event::Capture => [0, 0, 0],
                Event::MaxSteps => [0, 0, 0],
            };
            pix[0] = c[0];
            pix[1] = c[1];
            pix[2] = c[2];
        });

    let img = RgbImage::from_raw(width, height, buf)
        .context("failed to build image buffer")?;
    Ok(img)
}

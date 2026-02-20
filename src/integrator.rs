use crate::geodesic;
use crate::math;

#[derive(Debug, Clone, Copy)]
pub struct IntegratorParams {
    pub newton_iters: u32,
    pub newton_tol: f64,
}

fn f_mid(y_n: [f64; 8], z: [f64; 8], m: f64) -> [f64; 8] {
    let y_mid = math::scale8(math::add8(y_n, z), 0.5);
    geodesic::derivs(y_mid, m)
}

fn eval_f(z: [f64; 8], y_n: [f64; 8], h: f64, m: f64) -> [f64; 8] {
    let f = f_mid(y_n, z, m);
    let mut out = [0.0; 8];
    for i in 0..8 {
        out[i] = z[i] - y_n[i] - h * f[i];
    }
    out
}

pub fn implicit_midpoint_step(
    y_n: [f64; 8],
    h: f64,
    m: f64,
    params: IntegratorParams,
) -> Option<[f64; 8]> {
    let mut z = y_n;

    for _ in 0..params.newton_iters {
        let fz = eval_f(z, y_n, h, m);
        if !fz.iter().all(|v| v.is_finite()) {
            return None;
        }
        if math::inf_norm8(fz) < params.newton_tol {
            return Some(z);
        }

        let mut j = [[0.0; 8]; 8];
        for k in 0..8 {
            let eps = 1e-6 * (1.0 + z[k].abs());
            let mut zp = z;
            let mut zm = z;
            zp[k] += eps;
            zm[k] -= eps;
            let fp = eval_f(zp, y_n, h, m);
            let fm = eval_f(zm, y_n, h, m);
            let inv = 1.0 / (2.0 * eps);
            for i in 0..8 {
                j[i][k] = (fp[i] - fm[i]) * inv;
            }
        }

        let rhs = math::scale8(fz, -1.0);
        let dx = math::gauss_elim_solve(j, rhs)?;

        for i in 0..8 {
            z[i] += dx[i];
        }

        if !z.iter().all(|v| v.is_finite()) {
            return None;
        }
    }

    let fz = eval_f(z, y_n, h, m);
    if math::inf_norm8(fz) < params.newton_tol {
        Some(z)
    } else {
        None
    }
}

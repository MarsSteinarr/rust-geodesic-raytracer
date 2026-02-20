use crate::metric;

pub fn hamiltonian(x: [f64; 4], p: [f64; 4], m: f64) -> f64 {
    let g = metric::metric_inv(x, m);
    let mut s = 0.0;
    for a in 0..4 {
        for b in 0..4 {
            s += g[a][b] * p[a] * p[b];
        }
    }
    0.5 * s
}

pub fn derivs(y: [f64; 8], m: f64) -> [f64; 8] {
    let x = [y[0], y[1], y[2], y[3]];
    let p = [y[4], y[5], y[6], y[7]];

    let rr = metric::r(x);
    if rr < 1e-6 {
        return [f64::NAN; 8];
    }

    let g = metric::metric_inv(x, m);
    let dg = metric::metric_inv_derivs(x, m);

    let mut out = [0.0; 8];

    for mu in 0..4 {
        let mut v = 0.0;
        for nu in 0..4 {
            v += g[mu][nu] * p[nu];
        }
        out[mu] = v;
    }

    for mu in 0..4 {
        let mut s = 0.0;
        for a in 0..4 {
            for b in 0..4 {
                s += dg[mu][a][b] * p[a] * p[b];
            }
        }
        out[4 + mu] = -0.5 * s;
    }

    out
}

pub fn project_p0(x: [f64; 4], p_prev: [f64; 4], m: f64) -> Option<[f64; 4]> {
    let mut p = p_prev;

    let g = metric::metric_inv(x, m);
    let a = g[0][0];

    let mut b = 0.0;
    for i in 1..4 {
        b += 2.0 * g[0][i] * p[i];
    }

    let mut c = 0.0;
    for i in 1..4 {
        for j in 1..4 {
            c += g[i][j] * p[i] * p[j];
        }
    }

    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 || !disc.is_finite() {
        return None;
    }

    let sqrt_disc = disc.sqrt();
    if !sqrt_disc.is_finite() {
        return None;
    }

    let denom = 2.0 * a;
    if denom == 0.0 || !denom.is_finite() {
        return None;
    }

    let p0_1 = (-b + sqrt_disc) / denom;
    let p0_2 = (-b - sqrt_disc) / denom;

    let cand = [p0_1, p0_2];
    let mut best: Option<(f64, f64)> = None;

    for &p0 in &cand {
        let mut pp = p;
        pp[0] = p0;
        let mut dx0 = 0.0;
        for nu in 0..4 {
            dx0 += g[0][nu] * pp[nu];
        }
        if dx0 > 0.0 && dx0.is_finite() {
            let dist = (p0 - p_prev[0]).abs();
            match best {
                None => best = Some((dist, p0)),
                Some((best_dist, _)) if dist < best_dist => best = Some((dist, p0)),
                _ => {}
            }
        }
    }

    let (_, p0_best) = best?;
    p[0] = p0_best;
    Some(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn projection_enforces_null_constraint() {
        let m = 1.0;
        let x = [0.0, 0.3, -0.8, 5.0];
        let mut p = [-1.0, 0.2, 0.1, -0.3];
        p = project_p0(x, p, m).unwrap();
        let h = hamiltonian(x, p, m);
        assert!(h.abs() < 1e-10);
    }
}

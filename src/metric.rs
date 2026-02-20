pub fn r(x: [f64; 4]) -> f64 {
    (x[1] * x[1] + x[2] * x[2] + x[3] * x[3]).sqrt()
}

pub fn h_of_r(r: f64, m: f64) -> f64 {
    m / r
}

pub fn l_mu(x: [f64; 4]) -> [f64; 4] {
    let rr = r(x);
    if rr == 0.0 {
        [1.0, 0.0, 0.0, 0.0]
    } else {
        [1.0, x[1] / rr, x[2] / rr, x[3] / rr]
    }
}

pub fn l_up(x: [f64; 4]) -> [f64; 4] {
    let lm = l_mu(x);
    [-lm[0], lm[1], lm[2], lm[3]]
}

pub fn metric_inv(x: [f64; 4], m: f64) -> [[f64; 4]; 4] {
    let rr = r(x);
    let h = if rr == 0.0 { 0.0 } else { h_of_r(rr, m) };
    let lu = l_up(x);

    let mut g = [[0.0; 4]; 4];
    g[0][0] = -1.0;
    g[1][1] = 1.0;
    g[2][2] = 1.0;
    g[3][3] = 1.0;

    for a in 0..4 {
        for b in 0..4 {
            g[a][b] -= 2.0 * h * lu[a] * lu[b];
        }
    }

    g
}

pub fn metric_inv_derivs(x: [f64; 4], m: f64) -> [[[f64; 4]; 4]; 4] {
    let rr = r(x);
    let mut out = [[[0.0; 4]; 4]; 4];

    if rr < 1e-12 {
        return out;
    }

    let h = h_of_r(rr, m);
    let dhdx = {
        let r3 = rr * rr * rr;
        [0.0, -m * x[1] / r3, -m * x[2] / r3, -m * x[3] / r3]
    };

    let lu = l_up(x);

    let mut dl_up = [[0.0; 4]; 4];
    for sigma in 1..4 {
        for mu in 1..4 {
            let x_mu = x[mu];
            let x_sigma = x[sigma];
            let delta = if mu == sigma { 1.0 } else { 0.0 };
            dl_up[sigma][mu] = (delta - (x_mu * x_sigma) / (rr * rr)) / rr;
        }
    }

    for sigma in 0..4 {
        for a in 0..4 {
            for b in 0..4 {
                let term_h = -2.0 * dhdx[sigma] * lu[a] * lu[b];
                let term_l = -2.0 * h * (dl_up[sigma][a] * lu[b] + lu[a] * dl_up[sigma][b]);
                out[sigma][a][b] = term_h + term_l;
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn metric_inv_is_symmetric() {
        let x = [0.0, 1.3, -0.7, 2.1];
        let g = metric_inv(x, 1.0);
        for a in 0..4 {
            for b in 0..4 {
                let d = (g[a][b] - g[b][a]).abs();
                assert!(d < 1e-12);
            }
        }
    }

    #[test]
    fn metric_inv_dt_is_zero() {
        let x = [0.5, 2.0, 0.3, -1.0];
        let dg = metric_inv_derivs(x, 1.0);
        for a in 0..4 {
            for b in 0..4 {
                assert!(dg[0][a][b].abs() < 1e-12);
            }
        }
    }
}

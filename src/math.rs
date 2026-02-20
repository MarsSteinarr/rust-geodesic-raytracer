pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n == 0.0 || !n.is_finite() {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / n, v[1] / n, v[2] / n]
    }
}

pub fn add8(a: [f64; 8], b: [f64; 8]) -> [f64; 8] {
    let mut out = [0.0; 8];
    for i in 0..8 {
        out[i] = a[i] + b[i];
    }
    out
}

pub fn sub8(a: [f64; 8], b: [f64; 8]) -> [f64; 8] {
    let mut out = [0.0; 8];
    for i in 0..8 {
        out[i] = a[i] - b[i];
    }
    out
}

pub fn scale8(a: [f64; 8], s: f64) -> [f64; 8] {
    let mut out = [0.0; 8];
    for i in 0..8 {
        out[i] = a[i] * s;
    }
    out
}

pub fn inf_norm8(a: [f64; 8]) -> f64 {
    let mut m = 0.0;
    for i in 0..8 {
        let v = a[i].abs();
        if v > m {
            m = v;
        }
    }
    m
}

pub fn gauss_elim_solve(mut a: [[f64; 8]; 8], mut b: [f64; 8]) -> Option<[f64; 8]> {
    for k in 0..8 {
        let mut piv = k;
        let mut piv_abs = a[k][k].abs();
        for i in (k + 1)..8 {
            let v = a[i][k].abs();
            if v > piv_abs {
                piv_abs = v;
                piv = i;
            }
        }
        if piv_abs == 0.0 || !piv_abs.is_finite() {
            return None;
        }
        if piv != k {
            a.swap(k, piv);
            b.swap(k, piv);
        }

        let akk = a[k][k];
        for i in (k + 1)..8 {
            let f = a[i][k] / akk;
            a[i][k] = 0.0;
            for j in (k + 1)..8 {
                a[i][j] -= f * a[k][j];
            }
            b[i] -= f * b[k];
        }
    }

    let mut x = [0.0; 8];
    for i_rev in 0..8 {
        let i = 7 - i_rev;
        let mut s = b[i];
        for j in (i + 1)..8 {
            s -= a[i][j] * x[j];
        }
        let aii = a[i][i];
        if aii == 0.0 || !aii.is_finite() {
            return None;
        }
        x[i] = s / aii;
        if !x[i].is_finite() {
            return None;
        }
    }
    Some(x)
}

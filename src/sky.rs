use crate::math;

fn uv_from_dir(dir: [f64; 3]) -> Option<(f64, f64)> {
    let d = math::normalize3(dir);
    if d == [0.0, 0.0, 0.0] {
        return None;
    }
    let u = 0.5 + d[2].atan2(d[0]) / (2.0 * std::f64::consts::PI);
    let v = 0.5 - d[1].asin() / std::f64::consts::PI;
    Some((u, v))
}

pub fn sample_checker(dir: [f64; 3], squares: u32) -> [u8; 3] {
    let Some((u, v)) = uv_from_dir(dir) else {
        return [0, 0, 0];
    };

    let su = (u * squares as f64).floor() as i32;
    let sv = (v * squares as f64).floor() as i32;
    let parity = (su + sv) & 1;

    if parity == 0 {
        [255, 255, 255]
    } else {
        [40, 40, 40]
    }
}

fn hash_u32(mut x: u32) -> u32 {
    x ^= x >> 16;
    x = x.wrapping_mul(0x7feb352d);
    x ^= x >> 15;
    x = x.wrapping_mul(0x846ca68b);
    x ^= x >> 16;
    x
}

pub fn sample_stars(dir: [f64; 3]) -> [u8; 3] {
    let Some((u, v)) = uv_from_dir(dir) else {
        return [0, 0, 0];
    };

    let iu = (u * 4096.0).floor() as i32;
    let iv = (v * 2048.0).floor() as i32;

    let seed = ((iu as u32).wrapping_mul(0x9e3779b9)) ^ ((iv as u32).wrapping_mul(0x85ebca6b));
    let h = hash_u32(seed);

    let r = (h & 0xff) as u8;
    let g = ((h >> 8) & 0xff) as u8;
    let b = ((h >> 16) & 0xff) as u8;

    let base = [0u8, 0u8, 0u8];
    let p = (h >> 24) & 0xff;
    if p < 3 {
        let t = (p as f64) / 3.0;
        let intensity = (220.0 + 35.0 * (1.0 - t)) as u8;
        let warm = r as f64 / 255.0;
        let cool = b as f64 / 255.0;
        let rr = (intensity as f64 * (0.85 + 0.3 * warm)).clamp(0.0, 255.0) as u8;
        let gg = (intensity as f64 * (0.85 + 0.1 * g as f64 / 255.0)).clamp(0.0, 255.0) as u8;
        let bb = (intensity as f64 * (0.85 + 0.3 * cool)).clamp(0.0, 255.0) as u8;
        [rr, gg, bb]
    } else {
        base
    }
}

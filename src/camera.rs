use crate::math;

pub struct Camera {
    pub pos: [f64; 3],
    pub forward: [f64; 3],
    pub right: [f64; 3],
    pub up: [f64; 3],
    pub tan_half_fov: f64,
    pub aspect: f64,
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

pub fn make_camera(pos: [f64; 3], look_at: [f64; 3], up_hint: [f64; 3], fov_deg: f64, aspect: f64) -> Camera {
    let forward = math::normalize3([
        look_at[0] - pos[0],
        look_at[1] - pos[1],
        look_at[2] - pos[2],
    ]);
    let right = math::normalize3(cross(forward, up_hint));
    let up = math::normalize3(cross(right, forward));

    let fov_rad = fov_deg.to_radians();
    let tan_half_fov = (0.5 * fov_rad).tan();

    Camera {
        pos,
        forward,
        right,
        up,
        tan_half_fov,
        aspect,
    }
}

pub fn ray_dir(cam: &Camera, px: u32, py: u32, width: u32, height: u32) -> [f64; 3] {
    let x = (px as f64 + 0.5) / (width as f64);
    let y = (py as f64 + 0.5) / (height as f64);

    let ndc_x = 2.0 * x - 1.0;
    let ndc_y = 1.0 - 2.0 * y;

    let sx = ndc_x * cam.aspect * cam.tan_half_fov;
    let sy = ndc_y * cam.tan_half_fov;

    let dir = [
        cam.forward[0] + sx * cam.right[0] + sy * cam.up[0],
        cam.forward[1] + sx * cam.right[1] + sy * cam.up[1],
        cam.forward[2] + sx * cam.right[2] + sy * cam.up[2],
    ];

    math::normalize3(dir)
}

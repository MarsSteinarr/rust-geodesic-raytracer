use clap::Parser;
use clap::ValueEnum;

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum SkyMode {
    Black,
    Checker,
    Stars,
}

#[derive(Parser, Debug, Clone)]
#[command(name = "rust-geodesic-raytracer")]
#[command(about = "Null geodesic ray tracer around a Schwarzschild black hole (Kerr-Schild Cartesian)")]
pub struct Args {
    #[arg(long, default_value_t = 800)]
    pub width: u32,

    #[arg(long, default_value_t = 800)]
    pub height: u32,

    #[arg(long, default_value_t = 1.0)]
    pub m: f64,

    #[arg(long, default_value_t = 40.0)]
    pub fov_deg: f64,

    #[arg(long, value_names = ["X", "Y", "Z"], num_args = 3, allow_hyphen_values = true, value_parser = clap::value_parser!(f64))]
    pub cam_pos: Option<Vec<f64>>,

    #[arg(long, value_names = ["X", "Y", "Z"], num_args = 3, allow_hyphen_values = true, value_parser = clap::value_parser!(f64))]
    pub cam_look: Option<Vec<f64>>,

    #[arg(long, value_names = ["X", "Y", "Z"], num_args = 3, allow_hyphen_values = true, value_parser = clap::value_parser!(f64))]
    pub cam_up: Option<Vec<f64>>,

    #[arg(long, default_value_t = 0.03)]
    pub step_h: f64,

    #[arg(long, default_value_t = 20000)]
    pub max_steps: u32,

    #[arg(long, default_value_t = 8)]
    pub newton_iters: u32,

    #[arg(long, default_value_t = 1e-10)]
    pub newton_tol: f64,

    #[arg(long, default_value_t = 200.0)]
    pub r_sky: f64,

    #[arg(long, default_value_t = 2.0 * (1.0 + 1e-3))]
    pub r_cap: f64,

    #[arg(long, default_value_t = 16)]
    pub checker_squares: u32,

    #[arg(long, value_enum, default_value_t = SkyMode::Black)]
    pub sky: SkyMode,

    #[arg(long = "no-disk", default_value_t = true, action = clap::ArgAction::SetFalse)]
    pub disk: bool,

    #[arg(long, default_value_t = 6.0)]
    pub disk_r_in: f64,

    #[arg(long, default_value_t = 30.0)]
    pub disk_r_out: f64,

    #[arg(long, default_value_t = 0.08)]
    pub disk_half_thickness: f64,

    #[arg(long, default_value = "output.png")]
    pub outfile: String,

    #[arg(long, value_names = ["PX", "PY"], num_args = 2, value_parser = clap::value_parser!(u32))]
    pub debug_ray: Option<Vec<u32>>,
}

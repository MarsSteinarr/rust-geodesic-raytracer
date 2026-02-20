mod camera;
mod cli;
mod geodesic;
mod integrator;
mod math;
mod metric;
mod render;
mod sky;

use anyhow::Context;
use clap::Parser;

fn main() -> anyhow::Result<()> {
    let args = cli::Args::parse();

    if let Some(v) = args.debug_ray.as_deref() {
        if v.len() != 2 {
            anyhow::bail!("--debug-ray expects exactly 2 integers");
        }
        let px = v[0];
        let py = v[1];
        let dbg = render::trace_ray_debug(px, py, &args)?;
        println!("steps: {}", dbg.steps);
        println!("event: {}", dbg.event);
        println!("max_abs_h: {:.3e}", dbg.max_abs_h);
        println!("final_r: {:.6}", dbg.final_r);
        return Ok(());
    }

    let img = render::render_image(&args)?;
    img.save(&args.outfile)
        .with_context(|| format!("failed to save {}", args.outfile))?;
    Ok(())
}

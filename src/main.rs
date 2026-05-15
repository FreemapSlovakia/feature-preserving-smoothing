use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::Parser;
use rayon::prelude::*;
use wbgeotiff::{Compression, GeoTiff, GeoTiffWriter, SampleFormat, WriteLayout};

/// Reduce short-scale variation in a DEM using a modified
/// Sun et al. (2007) feature-preserving smoothing algorithm.
///
/// Port of the WhiteboxTools `FeaturePreservingSmoothing` tool.
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input raster DEM file (single-band GeoTIFF).
    #[arg(short = 'i', long = "dem")]
    input: PathBuf,

    /// Output raster file (GeoTIFF, float32).
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Filter kernel size (values < 3 are clamped to 3).
    #[arg(long, default_value_t = 11)]
    filter: usize,

    /// Maximum difference in surface normal vectors, in degrees.
    #[arg(long = "norm_diff", default_value_t = 15.0)]
    norm_diff: f32,

    /// Number of elevation-update iterations.
    #[arg(long = "num_iter", default_value_t = 3)]
    num_iter: usize,

    /// Maximum allowed absolute elevation change per cell.
    #[arg(long = "max_diff", default_value_t = 0.5)]
    max_diff: f32,

    /// Z conversion factor (for when vertical and horizontal units differ).
    /// If omitted, auto-derived for geographic CRSes, otherwise 1.0.
    #[arg(long = "zfactor")]
    zfactor: Option<f32>,
}

#[derive(Clone, Copy, Debug)]
struct Normal {
    a: f32,
    b: f32,
}

impl Normal {
    const ZERO: Normal = Normal { a: 0.0, b: 0.0 };

    /// Cosine of the angle between two surface normals (with implicit c = 1).
    #[inline]
    fn cos_angle(self, other: Normal) -> f32 {
        let denom = ((self.a * self.a + self.b * self.b + 1.0)
            * (other.a * other.a + other.b * other.b + 1.0))
            .sqrt();
        (self.a * other.a + self.b * other.b + 1.0) / denom
    }
}

fn main() -> Result<()> {
    let args = Args::parse();

    let filter_size = args.filter.max(3);
    let num_iter = args.num_iter.max(1);
    let norm_diff = args.norm_diff.clamp(0.0, 180.0);
    let max_z_diff = if args.max_diff.is_finite() {
        args.max_diff
    } else {
        f32::INFINITY
    };
    let threshold = norm_diff.to_radians().cos();

    println!("Reading {}", args.input.display());
    let start = Instant::now();

    let tiff = GeoTiff::open(&args.input)
        .with_context(|| format!("failed to open {}", args.input.display()))?;

    if tiff.band_count() < 1 {
        bail!("input has no bands");
    }
    let width = tiff.width() as usize;
    let height = tiff.height() as usize;
    if width < 3 || height < 3 {
        bail!("input is too small ({}x{}); need at least 3x3", width, height);
    }
    let nodata_in: Option<f32> = tiff.no_data().map(|v| v as f32);
    let nodata: f32 = nodata_in.unwrap_or(-32768.0);
    let input: Vec<f32> = tiff.read_band_f32(0).context("failed to decode band 0")?;
    if input.len() != width * height {
        bail!(
            "band 0 sample count {} does not match {}x{}",
            input.len(),
            width,
            height
        );
    }

    let geo_transform = tiff.geo_transform().cloned();
    let epsg = tiff.epsg();

    let (res_x, res_y, is_geographic) = match &geo_transform {
        Some(gt) => {
            let rx = gt.pixel_width.abs() as f32;
            let ry = gt.pixel_height.abs() as f32;
            let geo = epsg == Some(4326) || (rx <= 1.0 && ry <= 1.0 && epsg.is_none());
            (rx.max(f32::EPSILON), ry.max(f32::EPSILON), geo)
        }
        None => (1.0, 1.0, false),
    };

    let z_factor = match args.zfactor {
        Some(z) => z,
        None if is_geographic => {
            if let Some(gt) = &geo_transform {
                let mid_lat_deg = gt.y_origin + gt.pixel_height * (height as f64) / 2.0;
                let mid_lat = mid_lat_deg.clamp(-90.0, 90.0).to_radians();
                let z = (1.0 / (111_320.0 * mid_lat.cos())) as f32;
                println!("Input appears geographic; auto z-factor = {:.6e}", z);
                z
            } else {
                1.0
            }
        }
        None => 1.0,
    };

    let eight_res_x = res_x * 8.0;
    let eight_res_y = res_y * 8.0;

    let t_read = Instant::now();
    println!(
        "Read {}x{} band in {:.2}s",
        width,
        height,
        (t_read - start).as_secs_f64()
    );

    // 1. Surface normals (Horn 1981, 3x3).
    let normals = compute_normals(
        &input,
        width,
        height,
        nodata,
        z_factor,
        eight_res_x,
        eight_res_y,
    );
    let t_normals = Instant::now();
    println!(
        "Computed normals in {:.2}s",
        (t_normals - t_read).as_secs_f64()
    );

    // 2. Smooth the normal vector field.
    let smoothed = smooth_normals(&input, &normals, width, height, nodata, filter_size, threshold);
    let t_smooth = Instant::now();
    println!(
        "Smoothed normals in {:.2}s",
        (t_smooth - t_normals).as_secs_f64()
    );
    drop(normals);

    // 3. Iteratively update elevations.
    let output = update_elevations(
        &input,
        &smoothed,
        width,
        height,
        nodata,
        res_x,
        res_y,
        threshold,
        num_iter,
        max_z_diff,
    );
    let t_update = Instant::now();
    println!(
        "Updated elevations ({} iter) in {:.2}s",
        num_iter,
        (t_update - t_smooth).as_secs_f64()
    );

    println!("Writing {}", args.output.display());
    let mut writer = GeoTiffWriter::new(width as u32, height as u32, 1)
        .layout(WriteLayout::Tiled {
            tile_width: 256,
            tile_height: 256,
        })
        .compression(Compression::Deflate)
        .sample_format(SampleFormat::IeeeFloat)
        .software(format!(
            "feature-preserving-smoothing {}",
            env!("CARGO_PKG_VERSION")
        ));

    if let Some(gt) = geo_transform {
        writer = writer.geo_transform(gt);
    }
    if let Some(epsg) = epsg {
        writer = writer.epsg(epsg);
    }
    if let Some(nd) = nodata_in {
        writer = writer.no_data(nd as f64);
    }

    writer
        .write_f32(&args.output, &output)
        .with_context(|| format!("failed to write {}", args.output.display()))?;

    println!("Done in {:.2}s total", start.elapsed().as_secs_f64());
    Ok(())
}

/// Per-cell surface normals via Horn's (1981) 3x3 estimator. The vertical
/// component is implicitly 1, so only the x- and y-slopes are stored.
/// For neighbours that are out-of-bounds or nodata, the center's value is
/// substituted (matching the upstream WhiteboxTools behaviour).
fn compute_normals(
    input: &[f32],
    width: usize,
    height: usize,
    nodata: f32,
    z_factor: f32,
    eight_res_x: f32,
    eight_res_y: f32,
) -> Vec<Normal> {
    // NE, E, SE, S, SW, W, NW, N
    let dx: [isize; 8] = [1, 1, 1, 0, -1, -1, -1, 0];
    let dy: [isize; 8] = [-1, 0, 1, 1, 1, 0, -1, -1];

    let w = width as isize;
    let h = height as isize;

    let mut normals = vec![Normal::ZERO; width * height];

    normals
        .par_chunks_mut(width)
        .enumerate()
        .for_each(|(row, out_row)| {
            let row_i = row as isize;
            for col in 0..width {
                let col_i = col as isize;
                let z = input[row * width + col];
                if z == nodata {
                    continue;
                }
                let mut values = [0f32; 8];
                for i in 0..8 {
                    let yn = row_i + dy[i];
                    let xn = col_i + dx[i];
                    let zn = if yn >= 0 && yn < h && xn >= 0 && xn < w {
                        input[(yn as usize) * width + xn as usize]
                    } else {
                        nodata
                    };
                    values[i] = if zn != nodata { zn * z_factor } else { z * z_factor };
                }
                // p = [(z[NE] + 2 z[E] + z[SE]) - (z[NW] + 2 z[W] + z[SW])] / (8 dx)
                // q = [(z[NW] + 2 z[N] + z[NE]) - (z[SW] + 2 z[S] + z[SE])] / (8 dy)
                // Stored normal is (-p, -q, 1).
                let a = -(values[2] - values[4]
                    + 2.0 * (values[1] - values[5])
                    + values[0]
                    - values[6])
                    / eight_res_x;
                let b = -(values[6] - values[4]
                    + 2.0 * (values[7] - values[3])
                    + values[0]
                    - values[2])
                    / eight_res_y;
                out_row[col] = Normal { a, b };
            }
        });

    normals
}

/// Smooth the normal vector field with a bilateral-style filter: neighbours
/// whose normal is within `acos(threshold)` of the centre's are weighted by
/// `(cos_angle - threshold)^2`.
fn smooth_normals(
    input: &[f32],
    normals: &[Normal],
    width: usize,
    height: usize,
    nodata: f32,
    filter_size: usize,
    threshold: f32,
) -> Vec<Normal> {
    let midpoint = (filter_size / 2) as isize;
    let mut out = vec![Normal::ZERO; width * height];
    let w = width as isize;
    let h = height as isize;

    out.par_chunks_mut(width)
        .enumerate()
        .for_each(|(row, out_row)| {
            let row_i = row as isize;
            for col in 0..width {
                let col_i = col as isize;
                let z = input[row * width + col];
                if z == nodata {
                    continue;
                }
                let centre = normals[row * width + col];
                let mut sum_w = 0.0f32;
                let mut a = 0.0f32;
                let mut b = 0.0f32;
                for dr in -midpoint..=midpoint {
                    let yn = row_i + dr;
                    if yn < 0 || yn >= h {
                        continue;
                    }
                    let yi = yn as usize;
                    for dc in -midpoint..=midpoint {
                        let xn = col_i + dc;
                        if xn < 0 || xn >= w {
                            continue;
                        }
                        let xi = xn as usize;
                        let zn = input[yi * width + xi];
                        if zn == nodata {
                            continue;
                        }
                        let n = normals[yi * width + xi];
                        let cos_angle = centre.cos_angle(n);
                        if cos_angle > threshold {
                            let dw = cos_angle - threshold;
                            let weight = dw * dw;
                            sum_w += weight;
                            a += n.a * weight;
                            b += n.b * weight;
                        }
                    }
                }
                if sum_w > 0.0 {
                    out_row[col] = Normal { a: a / sum_w, b: b / sum_w };
                }
            }
        });

    out
}

/// Iteratively update elevations using the smoothed normal vector field.
///
/// At each cell, contributions are accumulated from 8 neighbours whose
/// smoothed normal is within `acos(threshold)` of the centre's. The update
/// is `z_new = sum(-(n.a*dx + n.b*dy - z_n) * w) / sum(w)`. Cells whose update
/// exceeds `max_z_diff` from the original elevation are reverted.
///
/// Mirrors the upstream Gauss-Seidel order: row by row, column by column,
/// reading already-updated values within the same iteration.
fn update_elevations(
    input: &[f32],
    smoothed: &[Normal],
    width: usize,
    height: usize,
    nodata: f32,
    res_x: f32,
    res_y: f32,
    threshold: f32,
    num_iter: usize,
    max_z_diff: f32,
) -> Vec<f32> {
    let dx: [isize; 8] = [1, 1, 1, 0, -1, -1, -1, 0];
    let dy: [isize; 8] = [-1, 0, 1, 1, 1, 0, -1, -1];
    let x_off = [-res_x, -res_x, -res_x, 0.0, res_x, res_x, res_x, 0.0];
    let y_off = [-res_y, 0.0, res_y, res_y, res_y, 0.0, -res_y, -res_y];

    let w = width as isize;
    let h = height as isize;
    let mut output: Vec<f32> = input.to_vec();

    for iter in 0..num_iter {
        for row in 0..height {
            let row_i = row as isize;
            for col in 0..width {
                let col_i = col as isize;
                let z_in = input[row * width + col];
                if z_in == nodata {
                    continue;
                }
                let centre_normal = smoothed[row * width + col];
                let mut sum_w = 0.0f32;
                let mut z_acc = 0.0f32;
                for n in 0..8 {
                    let yn = row_i + dy[n];
                    let xn = col_i + dx[n];
                    if yn < 0 || yn >= h || xn < 0 || xn >= w {
                        continue;
                    }
                    let zn = output[(yn as usize) * width + xn as usize];
                    if zn == nodata {
                        continue;
                    }
                    let neighbour_normal = smoothed[(yn as usize) * width + xn as usize];
                    let cos_angle = centre_normal.cos_angle(neighbour_normal);
                    if cos_angle > threshold {
                        let dw = cos_angle - threshold;
                        let weight = dw * dw;
                        sum_w += weight;
                        let predicted =
                            neighbour_normal.a * x_off[n] + neighbour_normal.b * y_off[n] - zn;
                        z_acc += -predicted * weight;
                    }
                }
                let new_z = if sum_w > 0.0 {
                    let candidate = z_acc / sum_w;
                    if (candidate - z_in).abs() <= max_z_diff {
                        candidate
                    } else {
                        z_in
                    }
                } else {
                    z_in
                };
                output[row * width + col] = new_z;
            }
        }
        println!("  iter {}/{}", iter + 1, num_iter);
    }

    output
}

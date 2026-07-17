# feature-preserving-smoothing

A small command-line tool that smooths a raster DEM while preserving breaks in
slope, using a modified version of the Sun et al. (2007) feature-preserving
mesh denoising algorithm.

The algorithm is ported (and lightly adapted) from the
[`FeaturePreservingSmoothing`](https://github.com/jblindsay/whitebox-tools/blob/master/whitebox-tools-app/src/tools/terrain_analysis/feature_preserving_smoothing.rs)
tool in [WhiteboxTools](https://www.whiteboxgeo.com/) by Dr. John Lindsay et al.
I/O is handled by the [`wbgeotiff`](https://crates.io/crates/wbgeotiff) crate.

This Rust implementation was written by [Claude](https://www.anthropic.com/claude).

## Algorithm

1. Compute a surface-normal vector per pixel from a 3×3 neighbourhood using
   Horn's (1981) estimator.
2. Smooth the normal vector field with a bilateral-style filter: neighbours
   whose normal is within `norm_diff` degrees of the centre's are weighted by
   `(cos_angle - cos_threshold)^2`.
3. Iteratively update elevations from the smoothed normals. Cells whose update
   would exceed `max_diff` from the original elevation are reverted, which
   preserves sharp features.

## Build

```bash
cargo build --release
```

## Usage

```text
Usage: feature-preserving-smoothing [OPTIONS] --dem <INPUT> --output <OUTPUT>

Options:
  -i, --dem <INPUT>            Input raster DEM file (single-band GeoTIFF)
  -o, --output <OUTPUT>        Output raster file (GeoTIFF, float32)
      --filter <FILTER>        Filter kernel size (values < 3 are clamped to 3) [default: 11]
      --norm_diff <NORM_DIFF>  Maximum difference in surface normal vectors, in degrees [default: 15]
      --num_iter <NUM_ITER>    Number of elevation-update iterations [default: 3]
      --max_diff <MAX_DIFF>    Maximum allowed absolute elevation change per cell [default: 0.5]
      --zfactor <ZFACTOR>      Z conversion factor (auto-derived for geographic CRSes, else 1.0)
  -h, --help                   Print help
  -V, --version                Print version
```

### Example

```bash
feature-preserving-smoothing \
  --dem input_dem.tif \
  --output smoothed_dem.tif \
  --filter 11 \
  --norm_diff 15 \
  --num_iter 3 \
  --max_diff 0.5
```

To increase the level of smoothing, raise `--norm_diff` (more neighbours
contribute) and/or `--num_iter`. The filter size has comparatively little
effect on the result.

## Notes

- Input must be a single-band GeoTIFF. The output is always written as a
  tiled, Deflate-compressed `Float32` GeoTIFF preserving the input's
  geo-transform, EPSG code, and `nodata` value.
- **Input must not use a TIFF predictor** (tag 317). `wbgeotiff` does not
  implement predictors: it inflates the compressed stream and returns the bytes
  verbatim, so `PREDICTOR=2`/`3` data would decode as garbage (`±Inf` /
  `f32::MAX`) with no error raised. The tool detects this and refuses to run.
  GDAL writes `PREDICTOR=3` routinely for float rasters — it compresses to
  roughly 1.7x versus 1.1x without — so this is easy to hit. Re-encode first:

  ```bash
  gdal_translate -co COMPRESS=DEFLATE -co PREDICTOR=1 predictor3.tif stripped.tif
  ```
- The first two stages (normal computation, normal smoothing) are parallelised
  with [`rayon`](https://crates.io/crates/rayon). The elevation-update stage is
  intentionally sequential (Gauss-Seidel style) to match the upstream
  WhiteboxTools behaviour.
- If the input is in geographic coordinates (EPSG:4326 or sub-degree pixels)
  and no `--zfactor` is given, a latitude-dependent factor is auto-derived
  using the same `1 / (111320 · cos(lat))` formula as WhiteboxTools.

## References

- Lindsay JB, Francioni A, Cockburn JMH. 2019. LiDAR DEM smoothing and the
  preservation of drainage features. *Remote Sensing*, 11(16), 1926.
  doi:10.3390/rs11161926
- Sun X, Rosin P, Martin R, Langbein F. 2007. Fast and effective
  feature-preserving mesh denoising. *IEEE Transactions on Visualization &
  Computer Graphics*, 13(5), 925–938.
- Horn BKP. 1981. Hill shading and the reflectance map. *Proceedings of the
  IEEE*, 69(1), 14–47.

## License

MIT.

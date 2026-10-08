# bactoloops.py -- the loop-calling macro in Python

A Python implementation of the Bactoloops macro (`../macro/Bactoloops.ijm`) for anyone who would
rather not run ImageJ/Fiji. It follows the macro step for step, with the same tunables (same names
and defaults), and writes the same outputs. See `../docs/macro_details.md` for what each
step does and why; this file only covers running the Python version and how closely it matches.

## Install

```
conda env create -f environment.yml
conda activate bactoloops
```

(or `pip install numpy scipy pandas tifffile numba matplotlib`). The first run compiles the
watershed code with numba, which takes a few seconds; after that a 461 x 461 map takes about
15 seconds.

## Run

```
python bactoloops.py map.tif --recenter-bin 368 --bin-size 10000 --out results/ --id wt_96h
```

- `map.tif`: the ICE-normalized contact map as a square float TIFF, one pixel per bin (as made by
  `../tools/hic_circular_recenter_to_tiff.py`).
- `--recenter-bin`: the `--recenter-bin` used to build the TIFF, or 0 if it wasn't recentered
  (368 for our *E. coli* maps, 0 for the Dugar et al. *B. subtilis* maps).
- `--bin-size`: bp per pixel. The tunables set in pixels were calibrated on 10 kb maps; see
  "Running at another resolution" in the macro README.

Outputs in `--out`, each ending in `_<id>`: `loop_calls` (CSV, same columns as the macro's table),
`Loop_density` and `Loop_density_smoothed` (16-bit TIFF), `HiC_APA` (TIFF with the log2 values, and
a PNG to view), `integrated_loop_density` (CSV), `Loop_density_smoothed_loops` (PNG with every loop
circled, ridge-search loops marked with a green +, the ridge mask outlined in cyan; low-confidence
loops are removed unless `keepLowConfidence = True`),
and `log` (the macro's Log lines). Tunables are the variables at the top of `bactoloops.py`.

## How it reproduces the macro

Rotation invariance (2026-10-06): the macro's rotation-invariant local background is
`ij_ops.periodic_background` (periodic tiling, all 16 rolling-ball grid phases, symmetrized), the
diagonal smoothing wraps around the map, and the donut kernels are oriented along the short arc.
Checked against the macro on *E. coli* WT 96 h (unrotated and rotated by 1, 2, 3, 57 and 230 bins)
and *B. subtilis* WT (unrotated, 97 and 211 bins): the same loops and confidence in all nine, every
score within 5e-4 (the macro's 3-decimal output); and the Python results are themselves identical
under rotation (to 1e-7).


The ImageJ operations the macro depends on are reimplemented from the ImageJ 1.54p source
(`ij_ops.py`, `ij_watershed.py`) rather than replaced by library equivalents, because the
library versions give different answers:

| ImageJ command | Python | Notes |
|---|---|---|
| Subtract Background (rolling ball, default options) | `ij_ops.rolling_ball_background` | ImageJ's algorithm: 3x3 mean pre-smoothing, image shrinking for large radii, ball patch, interpolation back up |
| Median | `ij_ops.median` | ImageJ's circular kernel (radius 2 = 21 pixels) |
| Gaussian Blur | `ij_ops.gaussian` | ImageJ's truncated, edge-corrected kernel; 8/16-bit results rounded |
| Size... (bilinear) | `ij_ops.resize_bilinear_16bit` | including ImageJ's quarter-scale centre shift |
| Binary > Watershed | `ij_watershed.watershed` | ImageJ's EDM and MaximumFinder segmentation, ported line for line |
| Analyze Particles | `ij_ops.particles` | 8-connected, in ImageJ's raster order, pixel-centre centroids |

Validation against the macro (clean Fiji runs of 2026-09-28) on all nine datasets
(*E. coli* WT 96 h, WT exponential, delta-hns, delta-stpA, delta-fis, delta-ihf; *B. subtilis* WT
run 1, WT runs 1+2, delta-rok runs 1+2):

- the same loops at the same positions in every dataset (0 to 214 loops per dataset);
- `confidence` identical for every loop except one delta-fis loop that sits on the 99th-percentile
  boundary (0.996 in the macro, 0.989 in Python);
- `log2_obs_exp`, `oe_percentile` and `integrated_excess` within 0.0005 (the macro writes three
  decimals) apart from that loop;
- `centroid1`/`centroid2` within 0.02 bins and `mean_intensity` within 0.2% -- a few particles gain
  or lose one pixel of the 5x-upsampled image at their edge;
- the rolling ball on its own, and the watershed on its own, are pixel-identical to ImageJ; the
  loop-density maps agree in 99.99% of pixels, the rest by 1 gray level (rounding).

Since 2026-09-30 (particles are found on the padded map, and the shot-noise statistics are measured on
the map without its padding; see "Circular edges" in the macro README), the two agree exactly on 10 of
12 maps (loop positions and confidence). On *E. coli* WT 96 h one particle is split in two by the
watershed in Python but not in ImageJ (125 vs 124 loops, 3 positions differ), and on delta-fis 2
positions differ. The cause is the noise statistics: measured on the unpadded map, the histogram's
lower edge is a pixel whose rolling-ball value differs between ImageJ and NumPy in the last float32
digits, which moves the noise floor by ~0.0005 and changes the integer rounding of a few pixels.

Remaining differences come from floating-point order of operations (ImageJ works in 32-bit float in
places where NumPy accumulates in 64-bit), which can flip near-ties.

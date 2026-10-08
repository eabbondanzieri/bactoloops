# Bactoloops macro: details

An ImageJ/Fiji macro (`macro/Bactoloops.ijm`; developed as `HiC_PeakFindingAlgorithm_with_symmetrization.ijm`) that takes an
ICE-normalized bacterial Hi-C contact map, removes the diagonal ridge and background to leave
a map of DNA-loop density, calls individual loops on that map, and reports each loop's
location in both TIFF-pixel and real genome coordinates. It extends the method described in
our preprint (PMID: 41279469).

**Python version.** `python/bactoloops.py` implements the same method in Python (no ImageJ
needed) and reproduces the macro's loop calls exactly on all our datasets; see
`python/README.md`.

## Requirements

- Fiji or ImageJ, with the image you want to analyze already open -- exported from an
  `.iced.matrix` file, at a fixed bin size such as 10kb. Either a 16-bit or a 32-bit TIFF
  works; a 32-bit source (able to hold negative values, if your export can produce them) is
  preferred, since the macro itself works in 32-bit through background subtraction and ridge
  removal regardless of the source's own bit depth -- see "Design notes" below.

## Running it

1. Open your Hi-C image in Fiji.
2. Check the tunables at the top of the macro (see below).
3. Run the macro. The startup dialog asks for the image to start from, a folder to save
   outputs to (pre-filled with ImageJ's current default directory -- change or clear it as
   needed), an identifier string, the **recenter bin**, the **bin size** (bp per pixel, default 10000), and a "Skip automatic saving" checkbox
   (see "Automatic output saving" below) -- check that box, or leave the folder blank, to skip
   automatic saving. The recenter bin must match the `--recenter-bin` used to build this TIFF
   with `tools/hic_circular_recenter_to_tiff.py` (0 if it wasn't recentered): 368 for our *E. coli*
   maps (ter centered, ori in the corners), 0 for the Dugar et al. *B. subtilis* GEO matrix. It
   is pre-filled from `recenterBin` at the top of the macro, and it only changes the `genome_*`
   columns and the `Integrated_loop_density` x-axis, not which loops are called. The value used
   is printed to the Log. It then works on a copy of the image you picked -- the original is never modified.
4. It's safe to re-run: any windows left over from a previous run are closed automatically
   before processing starts.

**Running it from the command line (batch runs).** Run Fiji headless, one dataset at a time:
`fiji --allow-multiple --no-splash -- --headless -macro run.ijm` (with this launcher
`--headless` goes after the `--`). The startup dialog can't be shown headless, so a batch script
sets `sourceTitle`, `saveDir`, `imageIdentifier`, `recenterBin`, `binSizeBp` and `skipSaving` in its
place, and saves the Log with `File.saveString(getInfo("log"), path)`. Without `--headless`, a
Fiji window that is already open (an interactive session) can receive the macro through
ImageJ's single-instance listener, even with `--allow-multiple`, and two macros running in one
Fiji share windows with the same names; in our runs that silently changed results between
otherwise identical runs. (`-port0`, which ImageJ documents for this, is ignored by this launcher.)
Headless runs gave byte-identical output on repeat.

## Pipeline overview

The step numbers match the `// ---- N.` section headers in the macro.

0. **Setup** -- pick the starting image and save options; duplicate it as `HiC_map` and
   convert to 32-bit.
1. **Scale and recolor** -- scale intensities so every condition's ICE row/column sum equals
   `targetRowSum` (see "Cross-condition scaling" below), and replace the main diagonal with a
   3x3 median to tame self-contact outliers.
2. **Pad** -- enlarge the canvas and wrap the signal periodically into the margin (the
   bacterial chromosome is circular), so later filters don't see a hard edge.
3. **Background subtraction** -- rolling-ball background estimate subtracted, giving
   `HiC_map_RB`.
4. **Noise measurement and ridge removal** -- measure a robust shot-noise median/SD on
   `HiC_map_RB` (on the map itself, not the padding), threshold a smoothed copy to build a mask of the diagonal ridge, and
   multiply it out, giving `Loop_density`.
5. **Noise floor** -- subtract the shot-noise floor from `Loop_density`, then convert it to
   16-bit (clipping negatives to 0).
6. **Smoothing** -- Gaussian-smooth a copy, giving `Loop_density_smoothed` (cropped back to the
   original size for display; the copy loops are called on stays padded through step 7).
7. **Loop calling** -- threshold an upsampled copy of the padded, smoothed map, split touching
   particles with a watershed, and measure each particle, keeping those whose centroid lies in the
   map itself rather than in the padding, so a loop at the edge of the map is found whole (7); collapse mirror-image calls on
   either side of the diagonal into one row per loop (7b); build a log2(observed/expected) map,
   subtract its large-scale local background, move each loop to the enrichment peak within
   one bin, rate every loop's confidence against it and drop the low-confidence ones, search the ridge-mask region with a
   HiCCUPS-style donut test, and measure each loop's integrated excess contacts (7c); draw each
   loop on `Loop_density_smoothed_loops` (7d); and average log2(observed/expected) contact
   frequency of `HiC_map` around every loop into `HiC_APA` (7e).
8. **Integrated loop density** -- crop `Loop_density` back to the original size and plot its
   column average against genome position, giving `Integrated_loop_density`.
9. **Save** -- write the six deliverables listed under "Automatic output saving".

## What it produces

All left open at the end for inspection/saving:

- **HiC_map_RB** -- background-subtracted map (32-bit). Step 4 measures the shot noise on
  THIS image, before the ridge is removed from it -- see "Design notes" below
- **Loop_density** -- `HiC_map_RB` with the diagonal ridge removed and the shot-noise floor
  subtracted. 32-bit until step 5, 16-bit after that
- **Loop_density_smoothed** -- the loop-density map WITHOUT the ridge mask, Gaussian-smoothed and
  cropped back to the original image size, with the ridge-mask region outlined (an overlay, kept
  in the saved TIFF). The one to use for figures. Loops are called on a copy with the ridge masked
  out; the ridge itself saturates at `displayMax`, so raise the display range in Fiji to look
  inside it
- **Loop_density_smoothed_outlines** -- outlines of the particles found by the loop scan, for a
  visual check of how well the watershed split touching particles apart
- **Loop_density_smoothed_outlines_labeled** -- the same outlines, with each particle's index
  number burned in (matches the `outline_ids` Results column)
- **Loop_density_smoothed_loops** -- an RGB copy of Loop_density_smoothed, upsampled by
  `loopCircleUpsampleFactor` (default 3x) so the circles render smoothly instead of blocky
  when zoomed in, with a circle drawn around every resolved loop, one circle per loop. The
  ridge-mask region, where the ridge search applies, is outlined in `ridgeOutlineColor` (cyan).
  Loops found by the ridge search are marked with a plus sign in `ridgeCallColor` (green).
  Low-confidence loops are removed by default (`keepLowConfidence`); if they are kept, they
  are circled in the lighter `lowConfidenceCircleColor`
- **HiC_APA** -- aggregate peak analysis (32-bit): log2(observed/expected) contact frequency,
  averaged over a small square of `HiC_map` centered on every resolved loop, shown on a
  red/blue LUT with white at 0 (see "APA" under "Design notes" below). A real,
  well-registered set of loop calls produces a symmetric, focused red spot in the middle of
  this image; a lopsided or off-center average is a sign something upstream is systematically
  biased (see "Design notes" below for the bin-indexing convention this is a useful check on).
  Enlarged `apaUpsampleFactor` times (nearest-neighbor) so it's easy to see; saved as a PNG.
- **HiC_map_log2OE** -- log2(observed/expected) of `HiC_map`, 3x3-averaged, with a
  rolling-ball local background (radius `localBackgroundRadius`, 500 kb) subtracted, on the same
  red/blue LUT as `HiC_APA` (range -1..1). The map the confidence columns are read from
  (step 7c)
- **Results table** -- one row per resolved loop, including any found by the ridge search,
  each with a confidence rating (see column reference below)
- **Integrated_loop_density** -- a plot of the integrated loop density (ILD, arbitrary
  units): the column average of `Loop_density`, as a
  continuous, genome-position-labeled track, for comparison against ChIP-seq/RNA-seq data.

### Results table columns

| Column | Meaning |
|---|---|
| `bin1`, `bin2` | Resolved loop anchors, as 0-based pixel/bin indices into *this TIFF*, after position refinement (step 7c) |
| `bp1`, `bp2` | `bin1`/`bin2` in bp (`bin * binSizeBp`) |
| `genome_bin1`, `genome_bin2` | `bin1`/`bin2` mapped back to the real genome bin, undoing the circular recenter if `recenterBin != 0`. Identical to `bin1`/`bin2` if `recenterBin == 0`. |
| `genome_bp1`, `genome_bp2` | `genome_bin1`/`genome_bin2` in bp -- the numbers to compare against ChIP-seq/RNA-seq tracks or a genome annotation |
| `n_detections` | 1 for a genuine single-sided detection, 2+ for a resolved mirror pair, 0 for a loop found by the ridge search |
| `note` | Which of those cases this row is: `mirror (clean or rounding-boundary)`, `TRUE single-sided detection`, or `ridge region (donut test)` |
| `outline_ids` | Index number(s) (comma-separated) of the particle(s) in Loop_density_smoothed_outlines_labeled this row came from |
| `mean_area` | Particle area, averaged across however many detections were merged |
| `mean_intensity` | Mean gray value, averaged the same way |
| `loop_strength` | `mean_area * mean_intensity` -- a single combined measure of how large and how intense the detection was. `outline_ids`, `mean_area`, `mean_intensity` and `loop_strength` are blank (NaN) for ridge-search loops, which have no particle; use `integrated_excess` to compare them with other loops |
| `log2_obs_exp` | The loop pixel's value in `HiC_map_log2OE`: 3x3-averaged log2(observed/expected) above the local background. A fold enrichment, comparable across datasets and loop distances |
| `oe_percentile` | Fraction of all pixels at the same genomic distance with a lower `log2_obs_exp` -- a distance-matched background |
| `confidence` | `high` if `oe_percentile >= confidenceHighPercentile` (0.99), `medium` if `>= confidenceMediumPercentile` (0.95), otherwise `low` |
| `centroid1`, `centroid2` | The particle centroids before rounding and refinement, in TIFF bins (pixel *i* spans *i* to *i* + 1), averaged over a mirror pair; blank for ridge-search loops. For checking how a call was placed |
| `refine_shift1`, `refine_shift2` | How far position refinement moved `bin1`/`bin2` (-1, 0 or +1 bins); 0 for ridge-search loops, which are placed on their own maximum |
| `integrated_excess` | Contacts above the local background, summed over a 5 x 5-bin (50 x 50 kb) window centered on the loop, as a percentage of a typical bin's non-local contacts. Roughly proportional to how often the loop forms; the recommended measure for comparing loop strengths, including between datasets |
| `donut_ratio`, `lowerleft_ratio`, `horizontal_ratio`, `vertical_ratio` | Observed / locally expected against the four HiCCUPS-style backgrounds (see "Ridge search" under "Design notes"). Reported for every loop; only ridge-search loops must pass them |

## Automatic output saving

The startup dialog's "Save outputs to" field is pre-filled with ImageJ's current default
directory, but it's optional -- clear it, or check "Skip automatic saving" (which overrides
the folder field either way, so you don't have to also clear or remember what's in it),
and nothing is saved automatically; all the windows and the Results table are still left
open at the end for a manual save. Otherwise, leave the folder as-is or
point it wherever you like (and give an identifier string) and, once the run finishes, six
files are written there automatically, each named with the identifier pasted onto the end:

| File | Contents |
|---|---|
| `Loop_density_<id>.tif` | The unsmoothed loop-density map (`loopDensityTitle`), 16-bit, cropped to the original image size -- the image `Integrated_loop_density` is computed from |
| `Loop_density_smoothed_loops_<id>.tif` | The upsampled loop-circle overlay image (`loopsTitle`) |
| `Loop_density_smoothed_<id>.tif` | The smoothed loop-density map (`loopDensitySmoothedTitle`) |
| `HiC_APA_<id>.png` | The aggregate peak analysis image (`apaTitle`), enlarged `apaUpsampleFactor` times and rendered with the red/blue LUT. The PNG is for viewing only -- read actual log2(obs/exp) values from the 32-bit `HiC_APA` window left open in Fiji (the center pixel's value is also printed to the Log) |
| `loop_calls_<id>.csv` | The Results table |
| `integrated_loop_density_<id>.csv` | The `Integrated_loop_density` plot data (`genome_position_bp,ILD_au`), written directly from the same arrays used to build the plot, rather than read back out of the Plot window |

The identifier is just pasted onto each filename as-is, so keep it filesystem-safe (letters,
numbers, `_`/`-`) -- something like a sample or condition name works well. Re-running with the
same folder and identifier overwrites that run's six files; use a different identifier (or
folder) to keep multiple runs side by side.

`saveAs("Tiff", ...)` renames the window it saves to match the saved filename as a side effect
(so `loopsTitle` briefly becomes something like `Loop_density_smoothed_loops_<id>` after being
saved) -- step 9 renames each of the three saved TIFF windows straight back to its canonical title
right after saving it, specifically so the `closeIfOpen()` calls at the very start of the macro
(which only know those canonical titles) still find and clear them on the next run. Without
that rename-back, a run that saved outputs would leave stale, oddly-named windows behind that
the next run's startup cleanup wouldn't recognize -- this way startup cleanup works the same way
whether or not saving happened on a previous run. `HiC_APA` doesn't need this: its PNG is saved from a disposable RGB
copy that is closed right after saving, so the `HiC_APA` window itself is never renamed.

## Coordinate conventions

- **`bin1`/`bin2` are 0-based.** A pixel at TIFF row/column index `i` corresponds to
  `bp = i * binSizeBp`, not `(i-1) * binSizeBp`.
- **Genome vs. TIFF coordinates.** If the TIFF was built with `--recenter-bin`, its own pixel
  grid is a circularly-rotated view of the genome -- `bin1`/`bin2` are useful for clicking
  back into the image itself, while `genome_bin1`/`genome_bin2`/`genome_bp1`/`genome_bp2` are
  what you actually want to compare against outside data.
- **The `Integrated_loop_density` x-axis can run negative.** When `recenterBin != 0`, the segment
  of the track that comes from before genome position 0 (wherever that lands in the TIFF) is
  shown as negative bp rather than wrapping back up near the genome's far end -- this keeps
  the whole track one continuous, monotonically increasing axis instead of a sawtooth.

## Tunables reference

**Cross-condition scaling** -- `calibrateTargetRowSum`, `targetRowSum`, `scalePercentile`,
`targetGrayLevel`: puts different Hi-C conditions on the same absolute intensity scale (the
same idea as "counts per million" for RNA-seq). `targetRowSum` should be chosen once against
a reference dataset and reused for every condition after that. To (re)calibrate: set
`calibrateTargetRowSum = true`, run the macro on the reference image, and read the suggested
`targetRowSum` printed to the Log window. Paste that value in, set
`calibrateTargetRowSum` back to `false`, and leave both alone from then on.

**Image processing** -- `contrastSaturated`, `padding`, `medianRadius`, `rollingBallRadius`,
`ridgeSmoothSigma`, `minParticleSize`, `maskBlurSigma`, `finalBlurSigma`: standard
background-subtraction/smoothing knobs; defaults work well for a 10kb bacterial-genome-scale
map and rarely need changing.

**Ridge/loop mask threshold** -- `ridgeThresholdMultiplier`: step 4 builds the ridge/loop mask
by thresholding a smoothed copy of `HiC_map_RB` at `shotNoiseMedian + ridgeThresholdMultiplier
* shotNoiseSD` (same `shotNoiseMedian`/`shotNoiseSD` described just below), via
`setRobustThreshold()` -- a fixed numeric cutoff rather than an automatic (e.g. Li) threshold.
Default `0.5`, a calibrated working value for this data/pipeline -- validated to reproduce the
Li threshold's mask and downstream loop-calling behavior (see "Design notes" below). Treat this
as a starting point for a new dataset and adjust based on how the mask looks against
`Loop_density`/`Loop_density_smoothed_outlines` for your data.

**Shot-noise subtraction** -- `noiseFloorMultiplier`: step 5 subtracts `shotNoiseMedian +
noiseFloorMultiplier * shotNoiseSD` from `Loop_density`, where `shotNoiseMedian`/
`shotNoiseSD` are a robust median and robust-SD (IQR / 1.349) measured on `HiC_map_RB`
(background-subtracted, but the ridge is still in it, and it's still 32-bit) rather than the
`Loop_density` itself -- see "Design notes" below for why that source
matters, not just the bit depth. Default `0.5`, a calibrated working value for this
data/pipeline -- lower than the old 16-bit pipeline's defaults needed, since the robust-SD is
now measured on clean, unclipped 32-bit statistics. Treat this as a starting point for a new
dataset (see `robustThresholdMultiplier` just below for how the two interact) and adjust based
on how much background noise is still getting through to
`Loop_density_smoothed_outlines`/`Loop_density_smoothed_loops`.

**Loop scan** -- `minSmallParticleArea`, `maxSmallParticleArea`, `upsampleFactor`,
`robustThresholdMultiplier`: controls the size range and threshold used to call individual
loop particles on an upsampled copy of the image. The threshold reuses `shotNoiseSD` from
step 4 (rather than recomputing a robust-SD from this working image's own histogram, which by
now has been shot-noise-subtracted and 16-bit-clipped -- exactly the skewed distribution
step 4's 32-bit measurement exists to avoid) at `(robustThresholdMultiplier -
noiseFloorMultiplier) * shotNoiseSD` above this image's own (already-shifted) zero point --
the subtraction accounts for step 5 having already shifted every pixel down by
`noiseFloorMultiplier * shotNoiseSD`, so the two multipliers are measuring loop-vs-noise
confidence on the same scale, one added on top of the other. `robustThresholdMultiplier` must
stay bigger than `noiseFloorMultiplier`, or every remaining (already background-subtracted)
pixel would pass as a loop. Default `1.5`, paired with `noiseFloorMultiplier`'s `0.5` above;
both multipliers are data-dependent, so treat `0.5`/`1.5` as starting points for a new dataset
and adjust based on how `Loop_density_smoothed_outlines`/`Loop_density_smoothed_loops` and the Results
table look for your data.

**Local background and integrated excess** -- `localBackgroundRadius` (50 bins = 500 kb): the
rolling-ball radius subtracted from the log2(obs/exp) map before confidence is scored; `0` turns
the correction off. `integratedExcessHalf` (2): half-width, in bins, of the window
`integrated_excess` sums over. See "Local background" under "Design notes" below.

**Position refinement** -- `refineRadius` (1 bin): each main-pass loop is moved to the highest
pixel of the single-pixel, background-corrected log2(obs/exp) map within this many bins; `0`
keeps the particle positions. See "Position refinement" under "Design notes" below.

**Loop confidence and ridge search** -- `confidenceHighPercentile` (0.99),
`confidenceMediumPercentile` (0.95): the `oe_percentile` cutoffs for the `confidence` column.
`keepLowConfidence` (false): loops rated `low` are removed from the table and the images right
after they are scored, before the ridge search, so they neither appear in the output nor block
ridge candidates next to them. Set it to true to keep them (for example, to inspect what the
95th-percentile cutoff removes).
`minLoopDist` (3 bins): no loop closer than this to the diagonal. `ridgeMinPercentile` (0.95),
`donutPeakWidth` (2), `donutWindow` (5), `donutRatioAll` (1.5), `donutRatioOne` (1.75): the ridge
search inside the ridge mask. See "Loop confidence and the ridge mask" under "Design notes" below.

**Display** -- `displayMax`, `loopCircleColor`, `ridgeCallColor`, `ridgeOutlineColor`, `ridgeOutlineWidth`, `lowConfidenceCircleColor`, `loopCircleRadius`, `loopCircleStrokeWidth`,
`loopCircleUpsampleFactor`, `apaSize`, `apaUpsampleFactor`, `apaColorLimit`: cosmetic settings
for the output images. `apaSize` must be odd. `apaUpsampleFactor` (default 40) enlarges
`HiC_APA` with nearest-neighbor scaling, so each averaged pixel becomes a solid block rather
than being interpolated. `apaColorLimit` sets `HiC_APA`'s color range to
`-apaColorLimit..+apaColorLimit` in log2(obs/exp) units; the default `0` picks the largest
absolute value in that APA automatically, so set a fixed value (e.g. `1`) when comparing
conditions side by side. `loopCircleUpsampleFactor` controls how much Loop_density_smoothed_loops is upsampled
before the circles are drawn, so they render smoothly rather than blocky when zoomed in --
set it to 1 for no upsampling.

**Genome coordinate mapping** -- `recenterBin`, `binSizeBp`: `recenterBin` is only the default
for the startup dialog's "Recenter bin" field, which is what a run actually uses -- enter exactly
the `--recenter-bin` value used to build this TIFF with `tools/hic_circular_recenter_to_tiff.py` (0 if
you didn't use that option). Values outside 0..(image size in bins) stop the macro with a
message. `binSizeBp` is likewise only the default for the dialog's "Bin size" field (bp per
pixel); enter the `--bin-size` the TIFF was built with.

**Running at another resolution.** The bin size sets every bp and genome coordinate, but most
tunables are set in pixels/bins and were calibrated on 10 kb maps: `padding`, `medianRadius`,
`rollingBallRadius`, `ridgeSmoothSigma`, `minParticleSize`, `maskBlurSigma`, `finalBlurSigma`,
`min/maxSmallParticleArea`, `apaSize`, `localBackgroundRadius` (50 bins = 500 kb),
`integratedExcessHalf` (2 bins = 50 x 50 kb window), `refineRadius`, `minLoopDist`,
`donutPeakWidth`, `donutWindow`, and the 3-bin duplicate distance in the ridge search. At another bin
size, decide for each whether it should keep its pixel value (image-processing scales) or keep its
physical size (e.g. halve the bin counts at 20 kb). The intensity scaling (`targetRowSum`) also
depends on bin size and needs recalibrating. The macro prints a reminder when the bin size isn't
10000 bp.

## Design notes

### Bin indexing and coordinate conversions

- **`bin1`/`bin2` are 0-based pixel indices**, so `bp = pixelIndex * binSizeBp` -- not
  `(pixelIndex - 1) * binSizeBp`, which `tools/hic_circular_recenter_to_tiff.py`'s own 1-based
  bin numbering might otherwise suggest. Validated by cross-correlating the macro's
  `Integrated_loop_density` output against an independently-built ChIP-seq enrichment z-score track across
  several candidate bin shifts: correlation peaks sharply (r=0.85, more than double any
  1-bin-shifted alignment) exactly at `bp = pixelIndex * binSizeBp`.
- **Particle centroids convert to bins with `floor()`, not `round()`.** ImageJ reports a
  particle's centroid in the convention where pixel index `i` spans the continuous range
  `[i, i+1)` -- so a particle centered on pixel 245 comes back as `X ~= 245.5`. `floor()` is
  the correct conversion back to a bin index; `round()` agrees with it only when the
  fractional part is below 0.5, and otherwise bumps the bin up by 1, biasing every affected
  loop toward smaller row *and* smaller column indices (visible as an off-center HiC_APA
  average). Validated the same way as bin indexing above: across 136 matched loops,
  switching to `floor()` raised mean anchor z-score from 2.95 to 3.17 and the fraction of
  loops with z>2 from 86% to 95% (paired test on the loops that moved, p=0.016).
- **Drawing a loop circle needs `bin + 0.5`, not `bin`.** The opposite conversion from the
  point above: `makeOval` needs to be centered on a pixel's visual *center* (`i + 0.5`), not
  its array *index* (`i`), so `Loop_density_smoothed_loops`'s circles are centered at
  `(bin + 0.5) * loopCircleUpsampleFactor`.

### APA

`HiC_APA` is built from `HiC_map` -- the scaled, diagonal-smoothed contact map, before any
background subtraction or ridge removal -- rather than from `Loop_density`, so it can be
compared directly with APA/pileup plots from other loop callers. It follows the
coolpup.py/cooltools convention:

1. **Expected** -- the mean contact frequency at each genomic distance, using the *circular*
   distance `d = min(|i-j|, N-|i-j|)` (N = number of bins), since the bacterial chromosome is
   circular.
2. **Observed/expected** -- every pixel in each loop's `apaSize x apaSize` window is divided by
   the expected value at its own distance, and these ratios are averaged across loops.
3. **log2** -- the average is shown as log2, so `0` (white) means no more contact than expected
   at that distance, red is enrichment and blue is depletion.

Pixels `<= 0` (bins removed by ICE) are treated as missing and skipped in both the expected and
the average. Windows near the edge of the map read the periodically wrapped signal from step 2's
padding. The Log reports the center pixel's log2(obs/exp) as a single summary number.

### Circular edges

The map is padded periodically in step 2, and every step that looks at a neighbourhood works on the
padded map, so no step sees a hard edge where the linear map starts: background subtraction, the
ridge mask, smoothing, particle finding, the obs/exp map and its local background, refinement,
the donut test, integrated excess and APA. Particle finding (step 7) runs on the whole padded map
and keeps only particles whose centroid lies in the map itself; each loop near an edge is then
found whole, once, rather than cut in two at the edge (until 2026-09-30 particles were found after
cropping, which lost a few loops within ~10 bins of the map edge: none in *E. coli* WT 96 h, where
the edge is at oriC and there are no H-NS loops, but several near oriC in *B. subtilis*). The
shot-noise statistics are measured on the map itself, because the padding repeats parts of it.

The output does not depend on where the circular map starts (since 2026-10-06). Four steps needed
care beyond the padding:

- **Local background (step 7c).** ImageJ's rolling ball shrinks the image 4-fold for radii of 30-100
  pixels, on a grid anchored at the image corner, so a plain Subtract Background depends on the map's
  start. With `localBackgroundPhaseAverage = true` (default) the 500-kb background is computed on a
  periodic 3 x 3 tiling of the map (cropped to the padded frame plus 100 px), averaged over all 16
  (x, y) grid phases (`localBackgroundPhases` = 4 per axis), and averaged with its transpose so it is
  exactly symmetric. Diagonal phases alone are not enough, because the genome length is generally not
  a multiple of 4, so the wrapped copies sit at different phases in x and y. Cost: about 5 s per map.
- **Diagonal smoothing (step 1).** Each diagonal pixel's 3x3 median wraps around the circular map
  (previously the two corner pixels used a truncated neighbourhood, which shifted the expected value
  at distance 0).
- **Ties in the distance-matched null.** A loop meets itself in its own null; with an exactly
  symmetric background these ties no longer depend on which side of the diagonal is read.
- **Donut kernels.** `donutRatios()` orients each pixel along the short arc of the chromosome, so the
  lower-left kernel always lies between the loop and the diagonal (previously it faced away from the
  diagonal for loops whose short arc crosses the map edge); a pixel exactly half the genome away gets
  the mean of both orientations.

Test: *E. coli* WT 96 h rotated by 1, 2, 3, 57, 118, 230, 231, 333 and 400 bins and *B. subtilis*
WT by 1, 2, 3, 97, 150, 211 and 300 bins give identical loops, confidence, log2(obs/exp),
`oe_percentile`, `integrated_excess` and donut ratios; only the particle descriptors
(`mean_area`, `mean_intensity`, `loop_strength`) differ, by at most 0.1% (rounding at the map edge
of the 5x-upsampled particles). Before these changes, rotations changed 3-9 of ~125 *E. coli* loops
and 2-3 of ~37 *B. subtilis* loops, all at the 95th-percentile boundary. Relative to the previous
version the calls change only at that boundary: WT 96 h 124 -> 128 loops (123 shared),
*B. subtilis* WT 37 -> 38 (36 shared); integrated excess r >= 0.995.

### Loop confidence and the ridge mask

Step 7c rates every loop against a distance-matched background in `HiC_map_log2OE`: `HiC_map`
divided by the expected contact frequency at each circular genomic distance (the same expected
the APA uses), log2-transformed, 3x3-averaged and corrected for local background (next section).
Pixels with no data (`<= 0`) are set to log2 = 0, so they neither create nor hide a peak.
`oe_percentile` is where the loop pixel falls among all pixels at the same distance, so a 300 kb
loop and a 1 Mb loop are judged on the same footing.

By default (`keepLowConfidence = false`) only medium- and high-confidence loops (at or above
the 95th percentile) are kept; the method works best restricted to these loops. The ridge search already applies the same 95th-percentile floor
(`ridgeMinPercentile`), so every loop in the output meets it.

### Local background

Expected-by-distance removes the fall-off of contact frequency with genomic distance and
nothing else. Any other large-scale structure stays in the log2(obs/exp) map and lifts every
pixel it covers, loop or not. In *B. subtilis* that is the SMC-mediated alignment of the two
chromosome arms (an X across the map) and the oriC and terminus domains; in *E. coli* the
large-scale structure is much weaker. Left in, it made loops that happen to sit on the X look
more significant than loops elsewhere, and made uncalled anchor pairs in *B. subtilis* look more
enriched than they are.

Step 7c therefore subtracts a rolling ball of radius `localBackgroundRadius` (50 bins = 500 kb;
Process > Subtract Background, default options) from the log2(obs/exp) map. 500 kb is large
compared with a loop (a few bins) and small compared with the domains and the arm-alignment X.
The map is still padded and periodically wrapped at this point, so the edges of the circular
chromosome are handled like the middle. Main-pass loop calls are unaffected; the correction
changes `log2_obs_exp`, `oe_percentile` and `confidence`, and which candidates the ridge search
accepts. It reproduces a Python implementation of the same correction
(`Supplemental_figure_Ecoli_vs_Bsub/localbg.py` in the paper-analysis repository, eabbondanzieri/HNS-Deep-Stationary-Bridging-2026) to within 0.001 log2 units.

`integrated_excess` uses the same local background. For each pixel, excess = observed -
expected x 2^(background the rolling ball removed), i.e. contacts above what distance plus the
large-scale structure predict (0 where the pixel has no data). It is summed over a 5 x 5-bin
window around the loop and expressed as a percentage of a typical bin's non-local contacts (the
median, over bins, of each bin's contacts >= 2 bins away; printed in the Log). `log2_obs_exp`
is a fold enrichment at the loop's centre pixel, which exaggerates long-range loops (their
expected is small) and misses loops spread over several bins; `integrated_excess` measures the
contacts the loop actually accounts for, which is roughly proportional to how often it forms.
Use `confidence` to decide whether a loop is real, and `integrated_excess` to compare how strong
loops are.

### Position refinement

A loop's position comes from the centroid of its particle on `Loop_density_smoothed`, and that
centroid lands on the loop's actual enrichment peak (the highest single-pixel log2(obs/exp) in its
3 x 3 neighbourhood) only about 40% of the time. Background subtraction, ridge removal, Gaussian
smoothing and the particle's shape each move it by a fraction of a bin, which at 10 kb is enough
to put the call one bin off. The loops are real, just placed coarsely: moving each call to its
local single-pixel maximum raised the anchors' mean ChIP-seq z-score, while moving each to a random
neighbour lowered it. The bin conversion (`floor()`) is not the cause. ImageJ's bilinear
`Size...` does shift the upsampled image by about 0.14 of an original pixel (a peak centred on
pixel *i* has its centroid at *i* + 0.36, not *i* + 0.5), but correcting for that did not move
more calls onto their peaks.

Step 7c therefore moves each main-pass loop to the highest pixel within `refineRadius` bins
(default 1, i.e. its 3 x 3 neighbourhood) of the single-pixel log2(obs/exp) map with the same
local background removed as `HiC_map_log2OE`. `refine_shift1`/`refine_shift2` record the move.
The 3 x 3-averaged map is deliberately not used for this: loops sit on anchor stripes, and
averaging pulls the maximum along the stripe, away from the dot itself (moving loops to the
averaged maximum lowered the anchors' ChIP z-score). Moves that would bring a loop closer than
`minLoopDist` to the diagonal are not allowed, and loops that land on the same pixel are
merged (this has not happened in our data).

Because the refined pixel is the best of its neighbourhood, a plain single-pixel null would
overstate its significance. `oe_percentile` for a refined loop is therefore computed against a
null built the same way: for every pixel at the same genomic distance, the same refinement pick
is made and the 3 x 3-averaged value read at the picked pixel. This makes `confidence` stricter
than before (WT 96 h: 125 of 146 loops medium or high, previously 148 of 150).

Tested against ChIP-seq and against six published loop callers (`Comparisons/`), medium- and
high-confidence loops, before and after refinement:

| | *E. coli* WT 96 h | *B. subtilis* WT (merged) |
|---|---|---|
| Mean anchor ChIP z-score (H-NS / Rok) | 2.37 -> 2.89 | 2.81 -> 3.55 |
| Mean z of the 10 strongest loops | 2.18 -> 3.11 | 3.41 -> 4.12 |
| APA P2LL | 1.43 -> 1.61 | 1.82 -> 2.08 |
| Other callers' loops at exactly the same pixel | 38% -> 60% | 32% -> 50% |

Delta-hns 96 h still gives no loops. The APA gain is partly circular (the refinement picks the
pixel with the highest obs/exp); the ChIP-seq gain is the independent evidence.

Loops close to the diagonal. Step 4's ridge mask removes every pixel connected to the diagonal
ridge, which also removes genuine loops close to it. The clearest example is the ~200 kb loop
between Rok clusters 2 and 3 in *B. subtilis* (Dugar et al. 2022, PMID 35075232), which the main
pass does not call. In our maps the mask covers everything within about 150 kb of the diagonal
and fades out by about 300 kb.

Ridge search. Inside the ridge mask a pixel can be high for its genomic distance without being
a loop -- domains and their corners lift whole neighbourhoods -- so step 7c holds candidates
there to a stricter, local standard (the approach HiCCUPS uses genome-wide, Rao et al. 2014):

1. **Candidates**: every local maximum (3 x 3) of the single-pixel, background-corrected
   log2(obs/exp) map that the position refinement uses, inside the ridge mask, at least
   `minLoopDist` (3) bins from the diagonal and not within 3 bins of a main-pass loop.
2. **Distance-matched percentile**: at least `ridgeMinPercentile` (0.95) against the same null
   refined loops are scored against.
3. **Donut test**: the candidate's contacts are compared with four local backgrounds -- the
   donut around it, the lower-left block towards the diagonal (the one a domain corner lifts),
   and a horizontal and a vertical strip through it, each window `donutWindow` (5) bins from the
   candidate with the central `donutPeakWidth` (2)-bin square left out. For each background,
   locally expected = (contacts summed over the background / distance-expected summed over it) x
   distance-expected at the candidate; pixels with no data or < 2 bins from the diagonal are left
   out. Observed / locally expected must be at least `donutRatioAll` (1.5) for all four and at
   least `donutRatioOne` (1.75) for the donut or the lower-left (the HiCCUPS thresholds at 10 kb).
   A real loop stands out from all four; a domain corner lifts its lower-left too and fails.
4. Passing candidates are accepted strongest first, skipping any within 3 bins of one already
   accepted. They are labelled `ridge region (donut test)` in the `note` column and marked with a
   plus sign on `Loop_density_smoothed_loops`, where the ridge-mask region is outlined.

The four ratios (`donut_ratio`, `lowerleft_ratio`, `horizontal_ratio`, `vertical_ratio`) are
reported for every loop as a check, but only ridge-search loops must pass them. The ratios are
not a significance test: HiCCUPS and cooltools also require a Poisson false-discovery rate, which
needs raw counts, and our *E. coli* maps are ICE-balanced. The false-positive rate was checked
empirically instead: of 12,000-15,000 random ridge pixels per map, none passed in *E. coli* WT
96 h, delta-hns or *B. subtilis* delta-rok (0.2% in *B. subtilis* WT, where random pixels land on
the real terminus loops).

Results (runs of 2026-09-29):

| Dataset | Ridge local maxima | >= 95th pct | Pass the donut test |
|---|---|---|---|
| *B. subtilis* WT stationary, runs 1 + 2 merged | 875 | 37 | 6: Rok cluster 2-yonX, 2-3 and yonX-3 (the three terminus-supercluster contacts Dugar et al. report; donut ratios 2.3-3.3), Rok cluster 7-8 (4.12-4.17 Mb, 5 bins from the diagonal), and the ~1.9 Mb site with yonX and with Rok cluster 3 |
| *B. subtilis* WT stationary, run 1 | 868 | 25 | 6: the same terminus contacts and Rok cluster 7-8, Rok cluster 1 with a site 180 kb away, and a contact next to the yonX site |
| *B. subtilis* delta-rok, runs 1 + 2 merged | 917 | 33 | 0 |
| *E. coli* WT 96 h | 791 | 24 | 0 |
| *E. coli* delta-hns 96 h | 846 | 27 | 0 |
| *E. coli* delta-stpA / delta-fis / delta-ihf 96 h | 691 / 684 / 440 | 23 / 18 / 14 | 2 / 1 / 4 |
| *E. coli* WT exponential | 865 | 36 | 2 |

Some ridge calls in the single-replicate *E. coli* maps have low or negative H-NS ChIP z-scores
(2.90-2.93 Mb in delta-stpA and delta-fis; 4.12-4.15 Mb in WT exponential, with donut ratios of
4-6): genuine focal short-range contacts that are not H-NS loops. The search needs an uncapped
matrix: on a TIFF exported with a display ceiling (e.g. clipped at 0.01) the terminus peaks
themselves are saturated.

This replaced an earlier anchor-guided search (an earlier development version, not included here),
which tested only pairs of loop anchors the main pass had already found, against the
distance-matched percentile alone. It found the same *B. subtilis* terminus contacts, but it
could only connect sites that already loop elsewhere, and it also accepted pairs that were high
for their distance without being focal (in *E. coli* WT 96 h its one addition has a donut ratio
of 1.28; in delta-rok it added a weak pair across oriC).

Alternatives tested (see `Dugar et al data/loop_network_analysis/ridge_threshold_sweep.*`):

- **A separate obs/exp peak scan over a fixed 50-300 kb band.** It recovered the same Rok
  cluster 2-3 terminus loops and added nothing in *E. coli*, but it was a second, open-ended
  detection method over part of the map, with an arbitrary distance window and no local
  (donut) test. The current ridge search scans the whole ridge mask but requires the donut test.
- **Raising `ridgeThresholdMultiplier` and relying on `confidence` to reject false positives.**
  Integrated loop density's ChIP correlation was unaffected (*E. coli* r = 0.844 at 0.5,
  0.847 at 5), and `confidence` did reject most of the new domain-contact calls. But in
  *B. subtilis* any value >= 1 also admits a run of calls on the SMC arm-alignment stripe near
  oriC, which is genuinely enriched over distance-matched expected and so passes `confidence`,
  while recovering Rok cluster 2-3 needed a value of 5. In *B. subtilis* the ridge mask at 0.5 is
  suppressing both domain contacts and arm alignment, so it stays at 0.5.

### 32-bit processing and noise measurement

The working image stays 32-bit float (rather than converting to 16-bit at the start) through
background subtraction and ridge removal, so those steps can produce genuine negative values
around zero instead of clipping at 0 the way an unsigned 16-bit image would. That matters
because `shotNoiseMedian`/`shotNoiseSD` (the robust median and IQR-based SD `computeRobustStats()`
measures at the start of step 4, before the ridge mask is built) need an unclipped,
symmetric-around-zero background distribution to describe correctly -- a distribution piled up
at exactly 0 would skew both statistics. Every step from there on -- the ridge/loop mask
threshold and the step 7 loop threshold alike -- reuses this one measurement instead of
recomputing a robust-SD from a later, already-subtracted and/or 16-bit-clipped image.

A few related details follow from that:

- The measurement is taken from `HiC_map_RB` (background-subtracted, ridge still in it),
  not `Loop_density` -- the ridge-removal mask built right after it zeroes out a
  chunk of pixels, which would itself skew a median/IQR estimate the same way clipping does.
  `HiC_map_RB` is still 32-bit and unclipped, and the ridge is a small enough fraction of
  the image that a robust estimator isn't bothered by it, the same way it isn't bothered by
  the sparse real loop peaks.
- Converting 32-bit to 16-bit uses ImageJ's default "Scale When Converting" behavior
  everywhere except step 5's one conversion, which wraps it in
  `setOption("ScaleConversions", false)`/`true` instead. Scaling remaps the image's current
  display range onto `[0, 65535]`, which is the right behavior for a conversion that just
  needs to preserve an image's overall shape but wrong for one that needs a literal, per-pixel
  clip at zero (negative -> 0), which is what "background clips at/below 0" requires here.
- Step 4's ridge/loop mask is thresholded with a manually-computed value
  (`shotNoiseMedian + ridgeThresholdMultiplier * shotNoiseSD`, via `setRobustThreshold()`)
  rather than an automatic method like ImageJ's Li threshold. An automatic, histogram-based
  threshold would need a 16-bit conversion first -- those methods bin a 32-bit image's
  histogram much more coarsely than a 16-bit one, which isn't fine enough for a clean
  ridge/background split -- but a manually-specified value has no such requirement, so the
  disposable `rbSmoothedTitle` copy used for the mask decision stays 32-bit like everything
  else in this stretch of the pipeline. Validated against an independent anti-HNS/WGS ChIP-seq
  enrichment track (see "Alternatives considered" below for how that comparison is used
  elsewhere in this pipeline's development): at `ridgeThresholdMultiplier = 0.5`, the resulting
  density track correlates with the prior Li-based mask's track at r=0.999, and the ChIP
  z-score correlation and loop-anchor enrichment are statistically indistinguishable from the
  Li-threshold pipeline.

See the "Known ImageJ macro quirks" section below for the scoping and string-concatenation
quirks `computeRobustStats()` works around to return its measurement out of a function.

### Alternatives considered

A couple of pipeline variants were tried and are worth documenting so they aren't
re-attempted without reason -- both were tested against an independent anti-HNS/WGS ChIP-seq
enrichment track (genome-wide track correlation, plus loop-anchor enrichment relative to that
track) and made no measurable difference on real data:

- **Deferring the 16-bit conversion until after blurring**, instead of converting
  `Loop_density` to 16-bit immediately after step 5's shot-noise subtraction. The
  idea was that blurring the real, unclipped, symmetric-around-zero 32-bit noise (rather than
  noise already floored to 0 by the 16-bit conversion) might leave `Loop_density_smoothed`'s
  background closer to a true zero. In practice this didn't change the loop-anchor ChIP
  enrichment; the pipeline keeps the simpler single, earlier conversion.
- **Upsampling before blurring**, instead of blurring `Loop_density_smoothed` at native
  resolution and only upsampling the working copies (`particlesTitle`/`grayTitle`) in
  step 7. The idea was that blurring at native resolution evaluates the Gaussian kernel on a
  coarse, blocky grid, and blurring on a finer, upsampled grid instead (with `finalBlurSigma`
  scaled by `upsampleFactor`) might give watershed cleaner, more isotropic peaks. This also
  made no measurable difference, so `upsampleFactor` remains a plain upsampling factor with
  no blur-order implications.

## Known ImageJ macro quirks this code works around

- With batch mode off, a `create` operation (`imageCalculator(... create)`, `Duplicate...`)
  followed by `rename()` or another command that acts on the active window can act on the wrong
  window, because the new window isn't always active yet. The macro therefore runs in batch mode
  from step 1 to the end and displays every window only when it's done
  (`setBatchMode("exit and display")`); batch-mode runs gave pixel-identical outputs to clean
  non-batch runs and are faster. As a second safeguard, the per-pixel maps in step 7c (the
  refinement map and the excess map behind `integrated_excess`) are computed in macro arrays read
  from each image by explicit title, and the obs/exp result window is selected by name before it
  is renamed.
- `return <expression> >= <value>;` returns the left-hand value and silently drops the
  comparison, so functions that return a true/false test use an explicit `if (...) return true;`.
  Calling code that then applies `!` to a number stops with an error dialog, which in headless
  mode has nowhere to show and hangs the run.
- A Fiji that is already open can receive a command-line macro through ImageJ's single-instance
  listener, even with `--allow-multiple`; run batch jobs headless (see "Running it").

- `&&` does not short-circuit on every ImageJ installation -- a compound condition that
  relies on the first clause protecting an out-of-bounds index on the second can still
  evaluate the second clause. `resolveSymmetricLoops()` uses nested `if`/`break` instead.
- The ternary `?:` operator doesn't reliably parse inside a string concatenation on every
  installation -- built as a plain `if`/`else` into a labeled variable instead where needed.
- `getHistogram()`'s 3-argument form only works for 8-bit/RGB images; 16-bit/32-bit images
  need the 5-argument form with an explicit range. Even then, the "values" output array isn't
  reliably populated -- code that needs it (`computeScaleFactor()`, `computeRobustStats()`)
  bins over the image's own actual `[min, max]` instead and reconstructs each bin's real value
  from its index and the known bin width, sidestepping "values" entirely. This also makes
  those functions work the same way on any bit depth, including a 32-bit image with negative
  values, rather than assuming a fixed 0-65535 range.
- Rewriting the Results table in place requires `run("Clear Results")`, not a plain row
  delete -- a row delete doesn't fully reset ImageJ's underlying measurement/Analyzer state,
  which otherwise leaves stale values in columns the macro didn't explicitly set.
- `roiManager("Show All with labels")` attaches a live overlay to the active window rather
  than drawing into pixels; `Flatten` bakes that overlay into a new image but does not clear
  it from the source window, and also switches the active window to the one it just created.
- Toggling `setBatchMode()` on and off repeatedly around individual steps is unreliable on
  some installations (windows can be silently dropped instead of just hidden) -- this macro
  turns batch mode on only for its first, most expensive step, and leaves it off after that.
- Reducing an image's bit depth (32-bit -> 16-bit, or down to 8-bit) *rescales* pixel values
  by default -- mapping the image's current display range onto the new type's full range --
  rather than clipping each pixel's own value in place. The display range at any given moment
  is leftover state from earlier steps, not necessarily related to the data currently being
  converted, so this can silently shift and rescale exactly the values a step is relying on
  being left alone (see "Design notes" above for how this applies to the step 5
  conversion). `setOption("ScaleConversions", false)` disables it for a literal clip; set back
  to `true` afterward for any later conversion that does want the default scaling behavior
  (this macro's RGB Color conversion in step 7d, for one).
- A function's internal variable assignments are local to that function, full stop -- even
  when a same-named variable already exists at the macro's top level. Pre-declaring a
  variable outside any function before calling a function that assigns to it does **not**
  turn the function's internal assignment into an update of that outer variable; the function
  silently creates and uses its own separate local copy instead, and the outer variable is
  left exactly as it was. The only mechanism confirmed reliable here for getting a value out
  of a function is a single `return` statement; `computeRobustStats()` returns both numbers
  packed into one comma-delimited string (`"robustMedian,robustSD"`), which the caller splits
  with `split()` and converts with `parseFloat()`.
- `+` does not reliably fall back to string concatenation just because one side is a string --
  it depends on which side. `"label=" + someNumber` (string first) concatenates correctly, and
  this pattern is used throughout the macro's `print()` calls without issue. But
  `someNumber + "," + someOtherNumber` (a NUMBER first) evaluates as arithmetic instead: the
  interpreter tries to parse `","` as a number, fails, and produces the numeric value `NaN`
  rather than a string -- which then breaks `split()` downstream (a `"NaN"` string has no
  comma to split on, so indexing past element 0 throws an "Index out of range" error). Leading
  with an empty string -- `"" + robustMedian + "," + robustSD` -- forces concatenation mode
  from the very first `+`. The rule of thumb for this codebase: always start a
  string-building `+` chain with a string literal, never a bare number.

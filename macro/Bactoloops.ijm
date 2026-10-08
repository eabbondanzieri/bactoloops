// ============================================================
//  HiC Loop Finding Algorithm
// ============================================================
//  Prompts you to pick which open image to start from, a folder to save
//  outputs to, and an identifier string, then works on a copy of the image
//  you picked -- it's never modified. Leave the save folder blank, or check
//  "Skip automatic saving", to skip automatic saving (step 9) entirely;
//  everything still stays open either way.
//
//  Works in 32-bit float through background subtraction and ridge removal
//  (steps 1-4) so those steps can go negative instead of clipping at 0 --
//  see the comments on steps 4 and 5 for why. Loop_density is converted to
//  16-bit partway through step 5, and every image derived from it after that
//  is 16-bit too (except HiC_APA, a 32-bit average, and
//  Loop_density_smoothed_loops, an RGB image).
//
//  Produces (left open at the end for inspection/saving):
//    - HiC_map_RB             background-subtracted map (32-bit) -- step 4
//                              measures the shot noise on THIS image, before
//                              the ridge is removed from it
//    - Loop_density            HiC_map_RB with the diagonal ridge removed
//                              and the shot-noise floor subtracted (32-bit
//                              until step 5, 16-bit from there on)
//    - Loop_density_smoothed   the same map WITHOUT the ridge mask applied,
//                              Gaussian-smoothed and cropped back to the
//                              original image size, with the ridge-mask region
//                              outlined (an overlay) -- the version for figures.
//                              Loops are called on a masked copy of it.
//    - Loop_density_smoothed_outlines
//                              outlines of the particles found in step 7,
//                              for a visual check of the watershed split
//    - Loop_density_smoothed_outlines_labeled
//                              the same outlines, with each particle's index
//                              number ("outline_ids" below) burned in
//    - Loop_density_smoothed_loops
//                              a copy of Loop_density_smoothed, upsampled
//                              loopCircleUpsampleFactor times, with a circle
//                              around every resolved loop (step 7d) -- one
//                              circle per loop, at whichever of its two
//                              symmetric (bin1,bin2)/(bin2,bin1) positions
//                              falls above the diagonal; loops rated "low"
//                              confidence (step 7c) are circled in a lighter
//                              color, loops found by the ridge search get a
//                              plus sign, and the ridge-mask region is outlined
//    - HiC_APA                 aggregate peak analysis -- log2(observed/
//                              expected) of HiC_map, averaged over an
//                              apaSize x apaSize square centered on every
//                              resolved loop (step 7e), on a red/blue LUT
//                              centered at 0; a real loop signal shows up as
//                              a red spot in the middle of this average. Enlarged
//                              apaUpsampleFactor times (nearest-neighbor) for
//                              viewing; saved as a PNG in step 9
//    - HiC_map_log2OE         log2(observed/expected) of HiC_map, 3x3
//                              averaged, with a rolling-ball local background
//                              subtracted (step 7c) -- the map the confidence
//                              columns are read from
//    - Results table           one row per resolved loop, including any added
//                              by the donut-tested ridge search inside the
//                              ridge mask (step 7c), each with a confidence rating
//                              and an integrated_excess (loop strength).
//                              See the README for the full column reference
//                              and coordinate conventions.
//    - Integrated_loop_density the column average of Loop_density (step 8),
//                              plotted against real genome position in bp,
//                              for comparison against ChIP-seq / RNA-seq
//                              tracks
//
//  If a save folder was given in the startup dialog, step 9 also writes six
//  files there automatically, each named with the identifier string pasted
//  onto the end: Loop_density_<id>.tif, Loop_density_smoothed_loops_<id>.tif,
//  Loop_density_smoothed_<id>.tif, HiC_APA_<id>.png, loop_calls_<id>.csv
//  (the Results table), and integrated_loop_density_<id>.csv (the
//  Integrated_loop_density plot data).
//
//  Safe to re-run: any windows left over from a previous run are
//  closed automatically before processing starts.
// ============================================================

// ---- window titles used by this macro --------------------------
// Each variable is named after the window it holds: loopDensityTitle holds
// "Loop_density", loopDensitySmoothedTitle holds "Loop_density_smoothed",
// and every window derived from Loop_density_smoothed is built from
// loopDensitySmoothedTitle plus a suffix that its variable name echoes
// (particlesTitle -> "_particles", outlinesTitle -> "_outlines", ...).
mainTitle                  = "HiC_map";
backgroundTitle            = "HiC_map_background";
rbTitle                    = "HiC_map_RB";                     // background-subtracted ("rolling ball")
rbSmoothedTitle            = "HiC_map_RB_smoothed";            // disposable smoothed copy used to build the ridge mask
maskTitle                  = "Mask of " + rbSmoothedTitle;      // ridge mask (named by Analyze Particles)
ridgeMaskTitle             = rbTitle + "_ridge_mask";         // copy of the ridge mask kept for step 7c's ridge search
                                                               // and for outlining the ridge on the display images
ridgeOutlineTitle          = ridgeMaskTitle + "_outline";      // disposable: ridge mask for building the outline
loopDensityCallTitle       = "Loop_density_smoothed_masked";   // disposable: Loop_density_smoothed with the ridge
                                                               // masked out -- what loops are called on (step 7)
loopDensityTitle           = "Loop_density";                    // ridge removed + shot-noise floor subtracted
loopDensitySmoothedTitle   = "Loop_density_smoothed";           // Gaussian-smoothed Loop_density
particlesTitle             = loopDensitySmoothedTitle + "_particles"; // disposable upsampled/binary/watershed working copy
grayTitle                  = loopDensitySmoothedTitle + "_gray";      // disposable upsampled grayscale copy, used only so
                                                                      // Mean gray value is sampled from loopDensitySmoothedTitle, not the mask
outlinesTitle              = loopDensitySmoothedTitle + "_outlines";
outlinesLabeledTitle       = outlinesTitle + "_labeled";        // outlinesTitle with particle index numbers burned in
loopsTitle                 = loopDensitySmoothedTitle + "_loops";     // loopDensitySmoothedTitle with a circle around each resolved loop
loopsSourceTitle           = loopsTitle + "_src";               // disposable working copy used to build loopsTitle
expectedTitle              = mainTitle + "_expected";          // disposable expected-by-distance image (step 7c)
obsExpTitle                = mainTitle + "_log2OE";            // log2(observed/expected), 3x3 averaged, local
                                                               // background subtracted (step 7c)
obsExpRawTitle             = obsExpTitle + "_distance_only";   // disposable: the same map before the local-
                                                               // background correction (step 7c)
rawLogTitle                = obsExpTitle + "_1px";             // disposable: log2(obs/exp) before the 3x3 average (step 7c)
apaTitle                   = "HiC_APA";                         // aggregate peak analysis image from step 7e
integratedLoopDensityTitle = "Integrated_loop_density";         // plot window from step 8 (column average of Loop_density)

// ---- tunable parameters ------------------------------------
// Cross-condition intensity scaling. targetRowSum puts every condition's
// ICE row/column-sum constant (see computeRowSum() below) on the same fixed
// number, the same idea as "counts per million" for RNA-seq libraries of
// different depths. Choose it ONCE against a reference dataset and reuse it
// for every condition after that -- changing it between runs defeats the
// purpose. To (re)calibrate: set calibrateTargetRowSum=true, run the macro
// on the reference image, and read the suggested targetRowSum the Log
// window prints (chosen so the scalePercentile-th percentile lands at
// targetGrayLevel out of 65535). Paste that value in below, set
// calibrateTargetRowSum back to false, and leave both alone after that.
calibrateTargetRowSum = false;    // true = print a suggested targetRowSum for THIS image; don't apply it yet
targetRowSum          = 5000000;  // <- replace with the value printed during calibration
scalePercentile       = 0.999;    // used only while calibrating -- see computeScaleFactor() below
targetGrayLevel       = 55000;    // used only while calibrating -- headroom below 65535 so the
                                   // brighter tail can clip instead of compressing the rest of the
                                   // range down to a handful of gray levels
contrastSaturated    = 0.35;
padding              = 200;    // total padding added around the plot (100 px/side)
medianRadius         = 2;
rollingBallRadius    = 10;
ridgeSmoothSigma     = 3;
minParticleSize      = 2000;   // px^2, minimum size kept by the ridge/loop mask
maskBlurSigma        = 1;
finalBlurSigma       = 1;
ridgeThresholdMultiplier = 0.5; // step 4 thresholds the ridge/loop mask at (shot-noise median) + this many
                                // robust-SD's above it, reusing the same shotNoiseMedian/shotNoiseSD
                                // computeRobustStats() measures from rbTitle -- replaces a Li auto-threshold.
                                // 0.5 is a calibrated working value for this data/pipeline, validated to
                                // reproduce the Li threshold's mask/loop-calling behavior on real data. Retune
                                // if you change datasets or upstream steps meaningfully, based on how the mask
                                // (or Loop_density) looks against your image.
noiseFloorMultiplier = 0.5;    // step 5 subtracts (shot-noise median) + this many robust-SD's above it -- a modest
                                // cutoff meant to knock out ordinary background noise while leaving real signal
                                // alone. Measured on rbTitle (32-bit, background-subtracted, ridge still in it --
                                // see step 4 and computeRobustStats() below for why that's the right source).
                                // 0.5 is a calibrated working value for this data/pipeline, not a generic default --
                                // it's deliberately low, since the robust-SD is now measured on clean, unclipped
                                // 32-bit statistics and doesn't need as large a multiplier as the old 16-bit
                                // pipeline did. Re-tune if you change datasets or upstream steps meaningfully.
displayMax           = 3000;
minSmallParticleArea = 4;      // px^2 at the ORIGINAL image scale -- lower size limit for the loop scan
maxSmallParticleArea = 100;     // px^2 at the ORIGINAL image scale -- upper size limit for the loop scan
upsampleFactor       = 5;      // linear upsampling applied before the loop scan
robustThresholdMultiplier = 1.5; // loop threshold, expressed as this many robust-SD's above the SAME shot-noise
                                // baseline noiseFloorMultiplier is measured against -- must be bigger than
                                // noiseFloorMultiplier, since a real loop call needs to clear a higher bar than
                                // "not noise". See step 7 and setRobustThreshold() below. 1.5 is a calibrated
                                // working value for this data/pipeline (paired with noiseFloorMultiplier=0.5
                                // above) -- re-tune both together if you change datasets or upstream steps
                                // meaningfully, based on how Loop_density_smoothed_outlines,
                                // Loop_density_smoothed_loops and the Results table look for your data.
refineRadius               = 1;     // step 7c: each main-pass loop is moved to the highest pixel of the corrected
                                    // log2(obs/exp) map within this many bins (1 = its 3 x 3 neighbourhood); the
                                    // particle centroid is only accurate to about a bin. 0 = keep particle positions
localBackgroundRadius      = 50;    // step 7c: rolling-ball radius (bins; 50 = 500 kb) subtracted from the
localBackgroundPhaseAverage = true;   // average the 500-kb background over ImageJ's 4 shrink-grid phases (rotation invariance)
localBackgroundPhases      = 4;
                                    // log2(obs/exp) map before confidence is scored, removing large-scale
                                    // structure that distance-only expected leaves in (e.g. the SMC arm
                                    // alignment and oriC/ter domains in B. subtilis). 0 = no correction
integratedExcessHalf       = 2;     // step 7c: integrated_excess sums a (2 x this + 1)^2-bin window around
                                    // each loop (2 -> 5 x 5 bins = 50 x 50 kb)
confidenceHighPercentile   = 0.99;  // "confidence" column: high at or above this oe_percentile ...
confidenceMediumPercentile = 0.95;  // ... medium at or above this, low below it
keepLowConfidence          = false; // false: loops rated "low" are removed from the table (and the images,
                                    // APA and anchors) right after they are scored; true keeps them
minLoopDist                = 3;     // bins -- no loop closer than this to the diagonal (position refinement
                                    // may not move a loop closer, and the ridge search starts here)
ridgeMinPercentile         = 0.95;  // step 7c ridge search: candidates inside the ridge mask must be at least
                                    // this distance-matched percentile ...
donutPeakWidth             = 2;     // ... and pass a HiCCUPS-style donut test (Rao et al. 2014; 10 kb defaults):
donutWindow                = 5;     //     peak half-width 2 bins, window half-width 5 bins,
donutRatioAll              = 1.5;   //     obs / local expected >= 1.5 against donut, lower-left, horizontal and
donutRatioOne              = 1.75;  //     vertical backgrounds, and >= 1.75 against the donut or the lower-left
loopCircleColor       = "blue"; // change if it doesn't read well against your image's actual LUT
ridgeCallColor        = "green"; // color of the plus marker for loops found by the ridge search (step 7c)
ridgeOutlineColor     = "#00c8e0"; // outline of the ridge mask (where the ridge search applies) on the display images
ridgeOutlineWidth     = 1.5;    // px, at the original image scale
loopCircleRadius      = 6;      // px, at loopDensitySmoothedTitle's (original-image) scale
loopCircleStrokeWidth = 1;      // px, at loopDensitySmoothedTitle's (original-image) scale
lowConfidenceCircleColor = "#99b8ff"; // lighter circle color for main-pass loops rated "low" confidence
                                 // (step 7c); set to loopCircleColor to draw them like the rest
loopCircleUpsampleFactor = 3;   // Loop_density_smoothed_loops is built on a copy of
                                 // loopDensitySmoothedTitle upsampled by this factor first, so the circles are drawn
                                 // with finer sub-pixel precision and look smooth rather than
                                 // blocky when zoomed in; 1 = no upsampling
apaSize               = 9;      // px, must be odd -- side length of the square averaged around
                                 // each loop for the aggregate peak analysis (step 7e) below
apaUpsampleFactor     = 40;     // HiC_APA is enlarged by this factor (nearest-neighbor, so each
                                 // averaged pixel becomes a solid block and no values are
                                 // interpolated) to make it easy to see -- 9 x 40 = 360 px square
apaColorLimit         = 0;      // HiC_APA color range, in log2(obs/exp) units, is -apaColorLimit..
                                 // +apaColorLimit (white = 0). 0 = automatic (the largest |value|
                                 // in this APA). Set a fixed value (e.g. 1) to put every condition
                                 // on the same color scale

// Genome coordinate mapping -- if this TIFF was built with
// hic_circular_recenter_to_tiff.py's --recenter-bin option, its bin 1 is NOT
// the genome's actual bin 1 -- the whole matrix was circularly rotated so a
// chosen bin landed in the middle of the image instead of at the corners.
// bin1/bin2/bp1/bp2 are always reported in THIS TIFF's own (possibly
// recentered) coordinates -- useful for clicking back into the image itself
// -- and genome_bin1/genome_bin2/genome_bp1/genome_bp2 undo that rotation
// (see undoRecenter() below) to give the REAL genome location, the numbers
// to actually compare against ChIP-seq/RNA-seq tracks or an annotation.
recenterBin = 368;      // DEFAULT for the "Recenter bin" field in the startup dialog, which is what
                       // the run actually uses. Enter EXACTLY the --recenter-bin value used to build
                       // THIS TIFF, or 0 if you didn't use --recenter-bin (then genome_bin1/genome_bin2
                       // just come out equal to bin1/bin2). 368 = our E. coli maps (ter centered,
                       // ori in the corners); the Dugar et al. B. subtilis GEO matrix is 0.
binSizeBp   = 10000;  // DEFAULT for the "Bin size" field in the startup dialog, which is what the run
                       // actually uses. bp per bin -- must match the --bin-size that run was given, since that's
                       // what both bp columns are computed from

// ---- 0. Pick the starting image & switch to a working copy ------
imageTitles = getList("image.titles");
if (imageTitles.length == 0)
    exit("Open the HiC image before running this macro.");

defaultSaveDir = getDirectory("current"); // ImageJ's own current default directory (e.g. wherever
                                           // the last file was opened/saved from) -- just a starting
                                           // point in the dialog below, not a resolved final path
Dialog.create("HiC Loop Finding");
Dialog.addChoice("Starting image:", imageTitles, imageTitles[0]);
Dialog.addDirectory("Save outputs to:", defaultSaveDir);
Dialog.addString("Image identifier:", "");
Dialog.addNumber("Recenter bin:", recenterBin, 0, 6, "(0 = TIFF not recentered)");
Dialog.addNumber("Bin size:", binSizeBp, 0, 8, "bp per pixel (must match the TIFF)");
Dialog.addMessage("Recenter bin must match the --recenter-bin used to build this TIFF\n" +
                  "(hic_circular_recenter_to_tiff.py). It only affects the genome_* columns\n" +
                  "and the Integrated_loop_density x-axis, not which loops are called.");
Dialog.addCheckbox("Skip automatic saving", false);
Dialog.show();
// Retrieved in the same order the fields were added above -- Dialog.getString()
// is also what returns an addDirectory() field's value, not a separate getter, and
// Dialog.getCheckbox() has to come after every addChoice()/addDirectory()/addString()
// /addNumber() getter, regardless of where addCheckbox() itself sat among those fields.
sourceTitle     = Dialog.getChoice();
saveDir         = Dialog.getString();
imageIdentifier = Dialog.getString();
recenterBin     = round(Dialog.getNumber());
binSizeBp       = round(Dialog.getNumber());
skipSaving      = Dialog.getCheckbox();
if (skipSaving)
    saveDir = "";   // step 9 already treats "" as "skip automatic saving" -- reuse that,
                     // regardless of whatever was typed/left in the folder field above
if (saveDir != "" && !endsWith(saveDir, File.separator))
    saveDir = saveDir + File.separator;

// Close anything left over from a previous run so this is safe to re-run.
closeIfOpen(mainTitle);
closeIfOpen("Result of " + mainTitle);
closeIfOpen(backgroundTitle);
closeIfOpen(rbTitle);
closeIfOpen(rbSmoothedTitle);
closeIfOpen(maskTitle);
closeIfOpen(ridgeMaskTitle);
closeIfOpen("Result of " + rbTitle);
closeIfOpen(loopDensityTitle);
closeIfOpen(loopDensitySmoothedTitle);
closeIfOpen(loopDensityCallTitle);
closeIfOpen(ridgeMaskTitle);
closeIfOpen(ridgeOutlineTitle);
closeIfOpen(particlesTitle);
closeIfOpen(grayTitle);
closeIfOpen(outlinesTitle);
closeIfOpen(outlinesLabeledTitle);
closeIfOpen(loopsTitle);
closeIfOpen(loopsSourceTitle);
closeIfOpen(apaTitle);
closeIfOpen(expectedTitle);
closeIfOpen(obsExpTitle);
closeIfOpen(obsExpRawTitle);
closeIfOpen(rawLogTitle);
closeIfOpen("Drawing of " + particlesTitle);
closeIfOpen(integratedLoopDensityTitle);
closeIfOpen(mainTitle + "_yplus");
closeIfOpen(mainTitle + "_yminus");
closeIfOpen(mainTitle + "_xplus");
closeIfOpen(mainTitle + "_xminus");
roiManager("reset");   // clear any ROIs left in the Manager from a previous run

setBatchMode(true);   // no window redraws while processing: faster, and see step 2 for why it stays on

// Work on a copy from here on; sourceTitle itself is never altered.
selectWindow(sourceTitle);
run("Select None");
run("Duplicate...", "title=" + mainTitle);
run("32-bit");   // Stay in 32-bit float through background subtraction and ridge removal (steps 1-4), so those
                 // steps can produce genuine negative values around zero instead of clipping at 0 the way a
                 // 16-bit (unsigned) image would -- that's what lets step 5 measure the true, unbiased shape of
                 // the background noise. Converted back to 16-bit partway through step 5, once that measurement
                 // is done; steps 6 onward are unaffected either way. Works the same regardless of whether the
                 // source TIFF itself is 16-bit or 32-bit.

// ---- 1. Scale intensities across conditions & recolor for viewing ------
// scaleFactor is derived from THIS image's own ICE row/column-sum constant
// (computeRowSum() below), so every condition's rows/columns end up summing
// to the same targetRowSum -- see the scaling note above the tunables.
rowSum = computeRowSum(mainTitle);

if (calibrateTargetRowSum) {
    // Calibration run: figure out what targetRowSum would need to be for
    // THIS image, and print it instead of guessing blindly -- still scales
    // this run normally, so a calibration run's own output still looks right.
    scaleFactor        = computeScaleFactor(mainTitle, scalePercentile, targetGrayLevel);
    suggestedTargetSum = rowSum * scaleFactor;
    print("CALIBRATION: this image's row-sum = " + rowSum + ". To put the " +
          (scalePercentile * 100) + "th percentile at gray level " + targetGrayLevel +
          " HERE, set targetRowSum = " + suggestedTargetSum +
          " and calibrateTargetRowSum = false, then reuse that targetRowSum for every condition.");
} else {
    scaleFactor = targetRowSum / rowSum;
}
print("HiC Loop Finding: row-sum = " + rowSum + ", scaleFactor = " + scaleFactor);

run("Multiply...", "value=" + scaleFactor);
run("Enhance Contrast", "saturated=" + contrastSaturated);
run("Fire");

imageSize    = getHeight();   // unpadded plot dimension
if (isNaN(recenterBin) || recenterBin < 0 || recenterBin > imageSize)
    exit("Recenter bin must be 0 (not recentered) or between 1 and " + imageSize +
         " (this TIFF's size in bins); got " + recenterBin + ".");
if (isNaN(binSizeBp) || binSizeBp <= 0)
    exit("Bin size must be a positive number of bp; got " + binSizeBp + ".");
print("HiC Loop Finding: " + sourceTitle + ", " + imageSize + " bins, recenter bin = " + recenterBin +
      ", bin size = " + binSizeBp + " bp");
if (binSizeBp != 10000)
    print("NOTE: the tunables that are set in pixels/bins (smoothing, particle sizes, ridge mask, local-background " +
          "radius, integrated-excess window, donut test, minimum loop distance) were calibrated on 10 kb maps -- see the README " +
          "before relying on calls at " + binSizeBp + " bp.");
negImageSize = -imageSize;
paddedSize   = imageSize + padding;
padHalf      = padding / 2;     // map bin b sits at padded pixel b + padHalf

// Replace every pixel on the central diagonal (x == y) with the median of
// its 3x3 neighborhood. Hi-C self-contact along the diagonal is often an
// extreme outlier compared to its immediate surroundings, which can throw
// off the background subtraction / ridge-removal / thresholding steps
// that follow if left as-is.
smoothDiagonal(mainTitle, imageSize);

// ---- 2. Pad the canvas and wrap the signal around the edges -----
// Adds shifted copies of the image into the new padding so later
// blurring / background subtraction don't see a hard, dark edge.
run("Canvas Size...", "width=" + paddedSize + " height=" + paddedSize + " position=Center");
wrapAndAdd(mainTitle, imageSize, negImageSize, "y");
wrapAndAdd(mainTitle, imageSize, negImageSize, "x");

// Batch mode stays ON for the rest of the macro: with windows displayed, a
// command that acts on the "active" image can run before a newly created
// window has become active, so it acts on the wrong image. That silently
// changed results between otherwise identical runs. Every image is shown at
// the very end (setBatchMode("exit and display")).

// ---- 3. Estimate & subtract background ---------------------------
selectWindow(mainTitle);
run("Duplicate...", "title=" + backgroundTitle);
run("Median...", "radius=" + medianRadius);
run("Subtract Background...", "rolling=" + rollingBallRadius + " create");

imageCalculator("Subtract create", mainTitle, backgroundTitle);
selectWindow("Result of " + mainTitle);
rename(rbTitle);
close(backgroundTitle);
// ---- 4. Measure shot noise, then build a loop/ridge mask -------------------
// Measured from rbTitle, NOT loopDensityTitle -- loopDensityTitle has a
// big island of exact-zero pixels everywhere the mask built below zeroes out
// the ridge, which is its own kind of skew (a pile-up at exactly 0, the same
// failure mode 32-bit was meant to avoid in the first place, just from
// masking instead of clipping). rbTitle is background-subtracted (step 3)
// but still has the ridge in it, so it's still 32-bit and unclipped --
// its background should genuinely hover around zero, with the ridge and any
// real loop peaks standing out. The ridge only covers a small fraction of
// the image, so a robust (median/IQR-based) estimate is insensitive to it
// the same way it's insensitive to the sparse loop peaks -- both are
// outliers a robust estimator is specifically designed to ignore.
// computeRobustStats() gets shotNoiseMedian/shotNoiseSD from this
// distribution, handed back as a single "median,SD" string (see that
// function's own comment for why it's returned rather than stored in a
// global) and split into two stably-named variables here so both the ridge
// mask threshold just below and step 7's loop threshold can reuse them,
// instead of recomputing a robust-SD from a later image's own histogram --
// which by then has been shot-noise-subtracted and/or 16-bit-clipped, and
// would suffer exactly the same skew problem this step exists to avoid.
// Measured on the unpadded map only: the padding holds wrapped copies of parts of the map, so including
// it would weight those parts twice (and make the thresholds depend on where the map happens to start).
// (computed on a cropped duplicate rather than through a selection, so the histogram is the plain
// whole-image one the Python version reproduces exactly)
selectWindow(rbTitle);
makeRectangle(padHalf, padHalf, imageSize, imageSize);
run("Duplicate...", "title=" + rbTitle + "_core");
selectWindow(rbTitle);
run("Select None");
robustStatsResult = computeRobustStats(rbTitle + "_core");
close(rbTitle + "_core");
robustStatsParts  = split(robustStatsResult, ",");
shotNoiseMedian   = parseFloat(robustStatsParts[0]);
shotNoiseSD       = parseFloat(robustStatsParts[1]);
print("Shot-noise estimate (32-bit " + rbTitle + ", pre-ridge-removal) -- median=" +
      shotNoiseMedian + ", robustSD=" + shotNoiseSD);

selectWindow(rbTitle);
run("Duplicate...", "title=" + rbSmoothedTitle);
run("Median...", "radius=" + medianRadius);
run("Gaussian Blur...", "sigma=" + ridgeSmoothSigma);

// Thresholds this smoothed copy at ridgeThresholdMultiplier robust-SD's above
// the shot-noise baseline just measured above, via setRobustThreshold() (the
// same helper step 7 uses below) -- a fixed numeric cutoff, replacing a Li
// auto-threshold. Stays 32-bit throughout: Li's histogram-based
// auto-detection needed a finer, 16-bit histogram to split reliably, but a
// manually-specified threshold value has no such requirement.
setRobustThreshold(ridgeThresholdMultiplier, shotNoiseMedian, shotNoiseSD);
// "redirect" is a global ImageJ preference that persists across separate
// macro runs -- set explicitly so a stale value from a previous run can't
// point at an image that doesn't exist yet in this one.
run("Set Measurements...", "area redirect=None decimal=3");
run("Analyze Particles...", "size=" + minParticleSize + "-Infinity show=Masks display clear");
selectWindow(maskTitle);
run("Gaussian Blur...", "sigma=" + maskBlurSigma);
run("Max...", "value=1");
run("Multiply...", "value=256");
run("Invert");
run("Divide...", "value=256");

imageCalculator("Multiply create", rbTitle, maskTitle);
selectWindow("Result of " + rbTitle);
rename(loopDensityTitle);
close(rbSmoothedTitle);
selectWindow(maskTitle);                 // kept until step 7c: 0 inside the ridge, 1 outside
run("Duplicate...", "title=" + ridgeMaskTitle);
close(maskTitle);

// ---- 5. Subtract residual shot noise -------
// Reuses shotNoiseMedian/shotNoiseSD, already measured above in step 4 from
// rbTitle (background-subtracted but still with the ridge in it, so still
// 32-bit and unclipped) -- no need to remeasure here.
selectWindow(loopDensityTitle);
noiseOffset     = shotNoiseMedian + noiseFloorMultiplier * shotNoiseSD;
print("Step 5 shot-noise floor: median=" + shotNoiseMedian + ", robustSD=" + shotNoiseSD +
      ", noiseFloorMultiplier=" + noiseFloorMultiplier + " -> subtracting " + noiseOffset);
run("Subtract...", "value=" + noiseOffset);

// loopDensityTitle converts to 16-bit right here, immediately after the
// shot-noise subtraction above, in a single conversion shared by both this
// window and loopDensitySmoothedTitle below (which inherits it via
// duplication). Nothing downstream needs loopDensityTitle to stay 32-bit any
// longer than this -- the step 8 Integrated_loop_density plot reads this
// same 16-bit, floor-clipped version. See the README for the tradeoffs considered
// in choosing this conversion point.
//
// ImageJ's default "Scale When Converting" behavior would rescale a 32-bit -> 16-bit
// conversion using the image's CURRENT DISPLAY RANGE (whatever min/max it happens to have
// inherited from earlier steps) mapped onto [0, 65535], instead of literally clipping each
// pixel's own value -- so a near-zero background wouldn't clip to 0 like it's supposed to; it
// would get linearly remapped to some arbitrary, uniformly-offset gray level instead. Turning
// scaling off makes the conversion a literal clip (negative -> 0, >65535 -> 65535) at each
// pixel's own value, which is what "background clips at/below 0" actually requires. Restored
// right after, since a later step (the RGB Color conversion in step 7d) DOES want
// display-range scaling.
setOption("ScaleConversions", false);
run("16-bit");
setOption("ScaleConversions", true);

// ---- 6. Loop_density_smoothed: smoothed copy for figures and loop calling ---
// loopDensitySmoothedTitle is duplicated from loopDensityTitle, which is
// already 16-bit and floor-clipped at this point, and the Gaussian blur below
// runs directly on that 16-bit copy -- no further bit-depth conversion needed afterward, since
// the duplicate and the blur both stay 16-bit throughout.
selectWindow(loopDensityTitle);
run("Duplicate...", "title=" + loopDensityCallTitle);
selectWindow(loopDensityCallTitle);
run("Gaussian Blur...", "sigma=" + finalBlurSigma);
// kept PADDED (periodically wrapped) through step 7, so a loop at the edge of the map is found whole;
// cropped after step 7
setMinAndMax(0, displayMax);
run("Invert LUT");

// The display version shows the data everywhere: the same steps, but from rbTitle (before the
// ridge mask was applied), with the ridge-mask region outlined in step 7d.
selectWindow(rbTitle);
run("Select None");
run("Duplicate...", "title=" + loopDensitySmoothedTitle);
run("Subtract...", "value=" + noiseOffset);
setOption("ScaleConversions", false);
run("16-bit");
setOption("ScaleConversions", true);
run("Gaussian Blur...", "sigma=" + finalBlurSigma);
run("Canvas Size...", "width=" + imageSize + " height=" + imageSize + " position=Center");
setMinAndMax(0, displayMax);
run("Invert LUT");

// ---- 7. Threshold & find loops on an upsampled copy ---------
// Works on an upsampled duplicate of Loop_density_smoothed (nearby loops are
// easier to separate at higher resolution); Loop_density_smoothed itself is
// never touched, so it's still exactly right for printing.
//
// A calibration is set on the upsampled copies so every measurement (area,
// centroid position) is reported back in the ORIGINAL image's pixel scale --
// so minSmallParticleArea/maxSmallParticleArea are still in original-image
// px^2, and centroid positions line up with the original image.
//
// Thresholds this (upsampled, already shot-noise-subtracted) working image
// at robustThresholdSD * shotNoiseSD above ITS OWN zero point -- reusing
// shotNoiseSD, the robust-SD computeRobustStats() measured back in step 4 on
// the clean 32-bit rbTitle, rather than recomputing a robust-SD from this
// image's own histogram. By this point the image has already been
// shot-noise-subtracted and 16-bit-clipped, which piles up a large fraction
// of background pixels at exactly 0 -- a badly skewed distribution to
// measure spread from, which is exactly the problem step 4's 32-bit
// measurement exists to avoid.
//
// robustThresholdSD = robustThresholdMultiplier - noiseFloorMultiplier: step
// 5 already subtracted shotNoiseMedian + noiseFloorMultiplier * shotNoiseSD
// from every pixel, so a pixel that was originally (in step 5's clean,
// pre-subtraction scale) robustThresholdMultiplier robust-SD's above
// shotNoiseMedian now sits, on THIS image's already-shifted scale, at
// exactly (robustThresholdMultiplier - noiseFloorMultiplier) * shotNoiseSD
// -- the algebra is: threshold_shifted = (shotNoiseMedian +
// robustThresholdMultiplier*shotNoiseSD) - (shotNoiseMedian +
// noiseFloorMultiplier*shotNoiseSD). robustThresholdMultiplier must stay
// bigger than noiseFloorMultiplier or this goes negative and every
// remaining (already-positive) pixel would pass as a "loop". Converts to a
// binary mask, then runs a watershed split to pull apart touching particles
// before measuring.
//
// Mean gray value is measured via a "redirect" image: grayTitle is a second
// upsampled copy of loopDensitySmoothedTitle that's never thresholded, so Analyze Particles
// finds particle shapes from the binary/watershed image (particlesTitle)
// but reads the Mean column from grayTitle's real grayscale intensities
// instead of the 0/255 mask values.
selectWindow(loopDensityCallTitle);
run("Select None");
run("Duplicate...", "title=" + particlesTitle);
upsampledSize = paddedSize * upsampleFactor;   // the padded map: loops at the map edges are found whole
run("Size...", "width=" + upsampledSize + " height=" + upsampledSize + " interpolation=Bilinear");
run("Set Scale...", "distance=" + upsampleFactor + " known=1 unit=px");

selectWindow(loopDensityCallTitle);
run("Select None");
run("Duplicate...", "title=" + grayTitle);
run("Size...", "width=" + upsampledSize + " height=" + upsampledSize + " interpolation=Bilinear");
// Calibrate grayTitle the same way as particlesTitle -- with "redirect"
// active, Area is reported using the redirect image's calibration, not the
// analyzed image's.
run("Set Scale...", "distance=" + upsampleFactor + " known=1 unit=px");

selectWindow(particlesTitle);
run("Select None");
run("Set Measurements...", "area mean centroid display redirect=" + grayTitle + " decimal=3");
robustThresholdSD = robustThresholdMultiplier - noiseFloorMultiplier;
setRobustThreshold(robustThresholdSD, 0, shotNoiseSD);
run("Options...", "iterations=1 count=1 black");   // objects=white(255) on black(0), as Watershed expects
run("Convert to Mask");
run("Watershed");

// show=[Bare Outlines] draws just the particle shapes, with no numbering of
// its own -- outline_id numbers are burned in separately below, on their
// own copy of this image.
//
// "add" pushes each particle into the ROI Manager, in the same order as the
// Results rows it also produces -- ROI Manager slot i and Results row i are
// the same particle. Renamed to a plain "0", "1", "2", ... below so that
// number is exactly what resolveSymmetricLoops() records as this particle's
// outline_id, and exactly what gets burned into outlinesLabeledTitle next.
run("Analyze Particles...", "size=" + minSmallParticleArea + "-" + maxSmallParticleArea + " show=[Bare Outlines] display clear add");
selectWindow("Drawing of " + particlesTitle);
rename(outlinesTitle);
close(particlesTitle);

// The scan ran on the padded, periodically wrapped map, so a loop near the map edge is found whole, and
// each loop can also appear a second time as a copy in the padding. Keep the one particle whose centroid
// lies in the central (unpadded) copy, and move centroids to unpadded pixel coordinates. Rows and ROI
// Manager slots are in the same order, so both are deleted together, last row first.
nScanned = nResults;
for (i = nResults - 1; i >= 0; i--) {
    cx = getResult("X", i); cy = getResult("Y", i);
    if (cx < padHalf || cx >= padHalf + imageSize || cy < padHalf || cy >= padHalf + imageSize) {
        IJ.deleteRows(i, i);
        roiManager("select", i);
        roiManager("delete");
    }
}
for (i = 0; i < nResults; i++) {
    setResult("X", i, getResult("X", i) - padHalf);
    setResult("Y", i, getResult("Y", i) - padHalf);
}
updateResults();
print("Step 7: " + nScanned + " particles on the padded map, " + nResults + " with their centroid inside the map.");

nOutlineRois = roiManager("count");
for (i = 0; i < nOutlineRois; i++) {
    roiManager("select", i);
    roiManager("rename", "" + i);
}
roiManager("deselect");

// Burn those index numbers onto a SEPARATE copy of the outlines image, so a
// loop's outline_id(s) in the Results table can be checked visually without
// cluttering the bare version. roiManager("Show All with labels") attaches a
// live overlay to the active window rather than drawing into pixels, and
// Flatten bakes that overlay into the new labeled image but doesn't clear it
// from the source -- so the source is explicitly reselected and cleared
// (both "Show None" and "Remove Overlay", belt-and-suspenders) afterward.
selectWindow(outlinesTitle);
roiManager("Show All with labels");
run("Flatten");
rename(outlinesLabeledTitle);

selectWindow(outlinesTitle);
roiManager("Show None");
run("Remove Overlay");
// both outline images were drawn on the padded, upsampled map: crop them to the map itself
for (k = 0; k < 2; k++) {
    if (k == 0) selectWindow(outlinesTitle); else selectWindow(outlinesLabeledTitle);
    makeRectangle(padHalf * upsampleFactor, padHalf * upsampleFactor, imageSize * upsampleFactor, imageSize * upsampleFactor);
    run("Crop");
    run("Select None");
}
selectWindow(loopDensityCallTitle);
run("Canvas Size...", "width=" + imageSize + " height=" + imageSize + " position=Center");

// Clear the redirect setting now that grayTitle is about to close, so a
// future run doesn't inherit a stale redirect pointing at a closed image.
run("Set Measurements...", "area mean centroid display redirect=None decimal=3");
close(grayTitle);

// ---- 7b. Resolve mirrored/duplicate loop calls across the diagonal --------
// The particle scan above finds each loop independently on each side of the
// diagonal, so almost every real loop shows up twice in the Results table,
// once as (X,Y) and once as (Y,X). This collapses those into one row per
// loop -- see resolveSymmetricLoops() below for how the collapsing works.
// Also maps every resolved loop back to real genome coordinates (undoing
// the circular recenter, if recenterBin != 0) via undoRecenter() below.
//
// bp1/bp2/genome_bp1/genome_bp2 are computed as bin * binSizeBp -- bin1/bin2
// are 0-based pixel indices into this TIFF. See the README for how this
// convention was established.
resolveSymmetricLoops(imageSize, recenterBin, binSizeBp);

// ---- 7c. Loop confidence & ridge search, from observed/expected --
// Rates every loop against a distance-matched background, using the contact
// map expressed as log2(observed/expected) rather than Loop_density.
// log2_obs_exp is the loop pixel's 3x3-averaged log2(obs/exp) after the local-
// background correction below; oe_percentile is the fraction of all pixels at
// the same genomic distance with a lower value, so short- and long-range loops
// are judged on the same footing. confidence is "high" at or above
// confidenceHighPercentile, "medium" at or above confidenceMediumPercentile,
// and "low" below that.
//
// Local-background correction: expected-by-distance removes the fall-off of
// contacts with genomic distance but nothing else, so large-scale structure
// (in B. subtilis the SMC arm alignment -- an X across the map -- and the
// oriC/ter domains) stays in the log2(obs/exp) map and lifts every pixel it
// covers, loop or not. A rolling ball of radius localBackgroundRadius bins
// (Process > Subtract Background, default options) is subtracted from the
// map to remove it. The map is still padded and periodically wrapped here, so
// the ball treats the edges of the circular chromosome like the middle.
//
// integrated_excess: contacts above the local background, summed over a
// (2 x integratedExcessHalf + 1)^2-bin window around the loop, as a
// percentage of a typical bin's non-local contacts (the median over bins of
// every contact >= 2 bins away). Per pixel, excess = observed - expected x
// 2^(background removed by the rolling ball), so it is measured against the
// same local background as log2_obs_exp. Unlike the fold enrichment in
// log2_obs_exp, it credits loops spread over several bins and is roughly
// proportional to how often the loop forms, which makes it the better
// measure for comparing loop strengths within or between datasets.
//
// Expected: the mean of every HiC_map pixel at the same CIRCULAR genomic
// distance d = min(|i-j|, N-|i-j|) (N = imageSize bins), since the bacterial
// chromosome is circular. Pixels <= 0 are treated as missing (bins ICE
// filtered out). mainTitle is still PADDED and periodically wrapped, so the
// log2(obs/exp) image is built at the padded size too and its 3x3 average
// wraps correctly at the edges; bin b sits at padded pixel b + padding/2.
// padHalf (= padding / 2) is set in step 2
nDist   = floor(imageSize / 2) + 1;       // circular distances 0..floor(N/2)

selectWindow(mainTitle);
run("Select None");
expSum   = newArray(nDist);
expCount = newArray(nDist);
rowNonLocal = newArray(imageSize);        // each bin's contacts >= 2 bins away (integrated_excess scale)
rowHasData  = newArray(imageSize);
for (r = 0; r < imageSize; r++) {
    for (c = r; c < imageSize; c++) {
        v = getPixel(c + padHalf, r + padHalf);
        if (v > 0) {
            d = minOf(c - r, imageSize - (c - r));
            expSum[d]   = expSum[d] + v;
            expCount[d] = expCount[d] + 1;
            rowHasData[r] = 1;
            rowHasData[c] = 1;
            if (d >= 2) {
                rowNonLocal[r] = rowNonLocal[r] + v;
                rowNonLocal[c] = rowNonLocal[c] + v;
            }
        }
    }
}
nonLocalList = newArray(imageSize);
nWithData = 0;
for (r = 0; r < imageSize; r++) {
    if (rowHasData[r] == 1) {
        nonLocalList[nWithData] = rowNonLocal[r];
        nWithData++;
    }
}
nonLocalList = Array.trim(nonLocalList, nWithData);
Array.sort(nonLocalList);
if (nWithData % 2 == 1)
    nonLocalMedian = nonLocalList[(nWithData - 1) / 2];
else
    nonLocalMedian = (nonLocalList[nWithData / 2 - 1] + nonLocalList[nWithData / 2]) / 2;
expected = newArray(nDist);
for (d = 0; d < nDist; d++) {
    if (expCount[d] > 0)
        expected[d] = expSum[d] / expCount[d];
}

// Expected-value image at the padded size, then obs/exp = mainTitle / it.
// Missing pixels (obs <= 0, or no expected) are set to obs/exp = 1, i.e.
// log2 = 0, "no information", so they neither create nor hide a peak.
newImage(expectedTitle, "32-bit black", paddedSize, paddedSize, 1);
for (y = 0; y < paddedSize; y++) {
    bi = ((y - padHalf) % imageSize + imageSize) % imageSize;
    for (x = 0; x < paddedSize; x++) {
        bj = ((x - padHalf) % imageSize + imageSize) % imageSize;
        d = abs(bi - bj);
        d = minOf(d, imageSize - d);
        if (expected[d] > 0)
            setPixel(x, y, expected[d]);
        else
            setPixel(x, y, 1);
    }
}
imageCalculator("Divide create 32-bit", mainTitle, expectedTitle);
selectWindow("Result of " + mainTitle);      // by name: the new window isn't always active yet
rename(obsExpTitle);
selectWindow(obsExpTitle);
run("Macro...", "code=[if (!(v > 0) || v > 1e30) v = 1;]");
run("Log");                                   // natural log ...
run("Multiply...", "value=" + (1 / log(2)));  // ... converted to log2
run("Duplicate...", "title=" + rawLogTitle);  // single-pixel copy, for position refinement below
selectWindow(obsExpTitle);
run("Convolve...", "text1=[1 1 1\n1 1 1\n1 1 1] normalize");   // 3x3 average

// Local-background correction (see above), keeping the uncorrected map
// just long enough to recover what the rolling ball removed.
run("Duplicate...", "title=" + obsExpRawTitle);
selectWindow(obsExpTitle);
if (localBackgroundRadius > 0) {
    if (localBackgroundPhaseAverage) {
        // ImageJ's rolling ball shrinks the image (4x at radius 50) on a grid anchored at pixel (0, 0), so its
        // background depends on where the circular map starts. The background is therefore computed on a
        // periodic copy of the map, averaged over every grid phase (localBackgroundPhases^2 shifts in x and y;
        // diagonal shifts alone are not enough, because the map length is generally not a multiple of 4, so
        // the wrapped copies sit at different phases in x and y). A rotation of the map then only shifts the
        // background. The periodic copy extends bgMargin px beyond the padded frame so the ball's own edge
        // effects (< 20 px) stay outside it.
        bgMargin = 100;
        bgSize   = paddedSize + 2 * bgMargin;
        selectWindow(obsExpTitle);
        makeRectangle(padHalf, padHalf, imageSize, imageSize);
        run("Duplicate...", "title=bgTile");
        selectWindow(obsExpTitle);
        run("Select None");
        selectWindow("bgTile");
        run("Canvas Size...", "width=" + (3 * imageSize) + " height=" + (3 * imageSize) + " position=Center zero");
        wrapAndAdd("bgTile", imageSize, negImageSize, "y");
        wrapAndAdd("bgTile", imageSize, negImageSize, "x");
        makeRectangle(imageSize - padHalf - bgMargin, imageSize - padHalf - bgMargin, bgSize, bgSize);
        run("Crop");
        newImage("bgPhaseSum", "32-bit black", bgSize, bgSize, 1);
        for (phy = 0; phy < localBackgroundPhases; phy++) {
            for (phx = 0; phx < localBackgroundPhases; phx++) {
                selectWindow("bgTile");
                run("Duplicate...", "title=bgPhase");
                selectWindow("bgPhase");
                if (phx + phy > 0) run("Translate...", "x=" + phx + " y=" + phy + " interpolation=None");
                run("Subtract Background...", "rolling=" + localBackgroundRadius + " create");
                if (phx + phy > 0) run("Translate...", "x=" + (-phx) + " y=" + (-phy) + " interpolation=None");
                imageCalculator("Add 32-bit", "bgPhaseSum", "bgPhase");
                close("bgPhase");
            }
        }
        close("bgTile");
        selectWindow("bgPhaseSum");
        run("Divide...", "value=" + (localBackgroundPhases * localBackgroundPhases));
        makeRectangle(bgMargin, bgMargin, paddedSize, paddedSize);   // back to the padded frame
        run("Duplicate...", "title=bgPhaseMean");
        close("bgPhaseSum");
        // The map is symmetric, so its background is too, up to float rounding in the phase sum; make it exactly
        // symmetric (average with its transpose) so that exact ties, e.g. a loop meeting itself in its own
        // distance-matched null, are not decided by which side of the diagonal a pixel is read from.
        selectWindow("bgPhaseMean");
        run("Duplicate...", "title=bgPhaseMeanT");
        run("Rotate 90 Degrees Right");
        run("Flip Horizontally");
        imageCalculator("Add 32-bit", "bgPhaseMean", "bgPhaseMeanT");
        close("bgPhaseMeanT");
        selectWindow("bgPhaseMean");
        run("Divide...", "value=2");
        imageCalculator("Subtract 32-bit", obsExpTitle, "bgPhaseMean");
        close("bgPhaseMean");
    } else
        run("Subtract Background...", "rolling=" + localBackgroundRadius);
}
close(expectedTitle);

// The per-pixel maps below are computed in macro arrays (index row * imageSize
// + column, unpadded), not with imageCalculator/rename chains: with batch mode
// off, a create-then-rename sequence occasionally acts on the wrong window,
// which silently changed integrated_excess (and position refinement) between
// otherwise identical runs. Each source image is read once, by explicit title.
//   rawArr    log2(obs/exp), single pixel (rawLogTitle)
//   avgRawArr the same, 3x3-averaged (obsExpRawTitle)
//   avgArr    3x3-averaged and local-background corrected (obsExpTitle)
//   obsArr    observed contacts (mainTitle)
nPix = imageSize * imageSize;
rawArr    = readUnpadded(rawLogTitle, imageSize, padHalf);
avgRawArr = readUnpadded(obsExpRawTitle, imageSize, padHalf);
avgArr    = readUnpadded(obsExpTitle, imageSize, padHalf);
obsArr    = readUnpadded(mainTitle, imageSize, padHalf);
close(rawLogTitle);
close(obsExpRawTitle);
// refArr: single-pixel log2(obs/exp) with the same local background removed,
// rawArr - (avgRawArr - avgArr), for position refinement.
// excArr: excess over the local background, observed - expected x 2^(removed
// background), 0 wherever observed is missing, for integrated_excess.
refArr = newArray(nPix);
excArr = newArray(nPix);
for (r = 0; r < imageSize; r++) {
    for (c = 0; c < imageSize; c++) {
        i = r * imageSize + c;
        bgRemoved = avgRawArr[i] - avgArr[i];
        refArr[i] = rawArr[i] - bgRemoved;
        if (obsArr[i] > 0)
            excArr[i] = obsArr[i] - expected[circDist(r, c, imageSize)] * pow(2, bgRemoved);
    }
}
selectWindow(obsExpTitle);

// Position refinement. A particle's centroid on Loop_density_smoothed lands on
// the loop's single-pixel enrichment peak only about 40% of the time --
// background subtraction, ridge removal, smoothing and the particle's shape
// each move it by a fraction of a bin -- so each main-pass loop is moved to the
// highest pixel within refineRadius bins of refArr: log2(obs/exp) at
// single-pixel resolution with the local background removed
// (refine_shift1/refine_shift2 record the move, in bins). The 3x3-averaged map
// is deliberately NOT used here: loops sit on anchor stripes, and averaging
// pulls the maximum along the stripe, away from the dot itself (tested against
// ChIP-seq: moving loops to the single-pixel maximum raised the anchors' mean
// H-NS/Rok z-score, moving them to the 3x3-averaged maximum lowered it). Moves
// that would bring a loop closer than minLoopDist to the diagonal are
// not allowed, and loops that land on the same pixel are merged into one row.
nMoved = 0;
nMergedRefine = 0;
if (refineRadius > 0) {
    for (row = 0; row < nResults; row++) {
        b1 = round(getResult("bin1", row));
        b2 = round(getResult("bin2", row));
        pick = refinePick(b1, b2, refArr, imageSize, refineRadius, minLoopDist);
        s1 = floor(pick / 1000) - 500;
        s2 = (pick % 1000) - 500;
        n1 = ((b1 + s1) % imageSize + imageSize) % imageSize;
        n2 = ((b2 + s2) % imageSize + imageSize) % imageSize;
        if (s1 != 0 || s2 != 0)
            nMoved++;
        setResult("bin1", row, minOf(n1, n2));
        setResult("bin2", row, maxOf(n1, n2));
        setResult("refine_shift1", row, s1);
        setResult("refine_shift2", row, s2);
    }
    // Merge loops that landed on the same pixel (keep the first, pool outline_ids).
    for (row = nResults - 1; row > 0; row--) {
        for (k = 0; k < row; k++) {
            if (getResult("bin1", k) == getResult("bin1", row) && getResult("bin2", k) == getResult("bin2", row)) {
                setResult("outline_ids", k, getResultString("outline_ids", k) + "," + getResultString("outline_ids", row));
                setResult("n_detections", k, getResult("n_detections", k) + getResult("n_detections", row));
                IJ.deleteRows(row, row);
                nMergedRefine++;
                k = row;                                              // stop scanning this row
            }
        }
    }
    for (row = 0; row < nResults; row++) {                         // bp and genome columns follow the new bins
        b1 = getResult("bin1", row);
        b2 = getResult("bin2", row);
        setResult("bp1", row, b1 * binSizeBp);
        setResult("bp2", row, b2 * binSizeBp);
        g1 = undoRecenter(b1, recenterBin, imageSize);
        g2 = undoRecenter(b2, recenterBin, imageSize);
        setResult("genome_bin1", row, g1);
        setResult("genome_bin2", row, g2);
        setResult("genome_bp1", row, g1 * binSizeBp);
        setResult("genome_bp2", row, g2 * binSizeBp);
    }
    updateResults();
    print("Step 7c refinement: " + nMoved + " of " + (nResults + nMergedRefine) + " loops moved to the single-pixel log2(obs/exp) maximum within " +
          refineRadius + " bin(s); " + nMergedRefine + " merged after landing on the same pixel.");
}
selectWindow(obsExpTitle);

// Confidence for every loop. A refined loop's position was chosen as the best
// of its neighbourhood, so its null is built the same way: for every pixel at
// the same genomic distance, the same refinement pick is made and the
// 3x3-averaged value read at the picked pixel.
nHigh = 0; nMedium = 0; nLow = 0;
for (row = 0; row < nResults; row++) {
    b1 = round(getResult("bin1", row));
    b2 = round(getResult("bin2", row));
    val = getPixel(b2 + padHalf, b1 + padHalf);
    d = circDist(b1, b2, imageSize);
    nBelow = 0;
    for (k = 0; k < imageSize; k++) {
        if (refineRadius > 0) {
            kc = (k + d) % imageSize;
            pick = refinePick(k, kc, refArr, imageSize, refineRadius, minLoopDist);
            pr = ((k + floor(pick / 1000) - 500) % imageSize + imageSize) % imageSize;
            pc = ((kc + (pick % 1000) - 500) % imageSize + imageSize) % imageSize;
            ref = avgArr[pr * imageSize + pc];
        } else
            ref = getPixel(((k + d) % imageSize) + padHalf, k + padHalf);
        if (ref < val)
            nBelow++;
    }
    pct = nBelow / imageSize;
    setResult("log2_obs_exp", row, val);
    setResult("oe_percentile", row, pct);
    if (pct >= confidenceHighPercentile) {
        setResult("confidence", row, "high");
        nHigh++;
    } else if (pct >= confidenceMediumPercentile) {
        setResult("confidence", row, "medium");
        nMedium++;
    } else {
        setResult("confidence", row, "low");
        nLow++;
    }
}
updateResults();
print("Step 7c confidence: " + nHigh + " high, " + nMedium + " medium, " + nLow + " low.");
if (!keepLowConfidence && nLow > 0) {
    for (row = nResults - 1; row >= 0; row--)
        if (getResultString("confidence", row) == "low")
            IJ.deleteRows(row, row);
    updateResults();
    print("Step 7c: removed the " + nLow + " low-confidence loop(s) (keepLowConfidence = false); " + nResults + " remain.");
}

// Ridge search. The ridge mask (step 4) blanks everything connected to the
// diagonal ridge -- mostly short-range contacts and chromosomal domains, where a
// pixel can be high for its distance without being a loop -- so the main pass
// can't call loops there. Inside the mask, loops are therefore held to a stricter,
// local standard: every local maximum of refArr (the single-pixel,
// background-corrected map the refinement uses) at least minLoopDist bins from
// the diagonal and not within 3 bins of a main-pass loop is a candidate. It must
// reach ridgeMinPercentile of its distance-matched null (built as for refined
// loops) AND pass a HiCCUPS-style donut test (Rao et al. 2014): its contacts
// divided by the contacts expected from each of four local backgrounds -- the
// donut, the lower-left block towards the diagonal (which catches domain
// corners), a horizontal and a vertical strip, each scaled by the
// distance-expected -- must all be at least donutRatioAll, and at least
// donutRatioOne for the donut or the lower-left. Passing candidates are accepted
// strongest first, skipping any within 3 bins of one already accepted.
ridgeArr = readUnpadded(ridgeMaskTitle, imageSize, padHalf);
nullCache = newArray(imageSize * nDist);
nullDone  = newArray(nDist);
donutOut  = newArray(4);
nRidgeMax = 0; nRidgePct = 0; nRidgeCand = 0; nRidgeAdded = 0;
if (nResults > 0) {
    nMain  = nResults;
    mainB1 = newArray(nMain);
    mainB2 = newArray(nMain);
    for (i = 0; i < nMain; i++) {
        mainB1[i] = round(getResult("bin1", i));
        mainB2[i] = round(getResult("bin2", i));
    }
} else {
    nMain = 0; mainB1 = newArray(1); mainB2 = newArray(1);
}
candVal = newArray(imageSize * nDist);
candB1  = newArray(imageSize * nDist);
candB2  = newArray(imageSize * nDist);
candPct = newArray(imageSize * nDist);
for (i = 0; i < imageSize; i++) {
    for (d = minLoopDist; d < nDist; d++) {
        j = (i + d) % imageSize;
        if (ridgeArr[i * imageSize + j] >= 0.5)
            continue;                                   // outside the ridge: the main pass covers it
        v = refArr[i * imageSize + j];
        isMax = true;
        for (dy = -1; dy <= 1; dy++)
            for (dx = -1; dx <= 1; dx++)
                if (refArr[(((i + dy) % imageSize + imageSize) % imageSize) * imageSize + (((j + dx) % imageSize + imageSize) % imageSize)] > v)
                    isMax = false;
        if (!isMax)
            continue;
        nRidgeMax++;
        b1 = minOf(i, j);
        b2 = maxOf(i, j);
        near = false;
        for (k = 0; k < nMain; k++)
            if (abs(b1 - mainB1[k]) <= 3 && abs(b2 - mainB2[k]) <= 3)
                near = true;
        if (near)
            continue;                                   // already a main-pass loop
        if (nullDone[d] == 0) {                         // distance-matched null, built as for refined loops
            for (k = 0; k < imageSize; k++) {
                kc = (k + d) % imageSize;
                pick = refinePick(k, kc, refArr, imageSize, refineRadius, minLoopDist);
                pr = ((k + floor(pick / 1000) - 500) % imageSize + imageSize) % imageSize;
                pc = ((kc + (pick % 1000) - 500) % imageSize + imageSize) % imageSize;
                nullCache[d * imageSize + k] = avgArr[pr * imageSize + pc];
            }
            nullDone[d] = 1;
        }
        val = avgArr[i * imageSize + j];
        nBelow = 0;
        for (k = 0; k < imageSize; k++)
            if (nullCache[d * imageSize + k] < val)
                nBelow++;
        pct = nBelow / imageSize;
        if (pct < ridgeMinPercentile)
            continue;
        nRidgePct++;
        donutRatios(i, j, obsArr, expected, imageSize, donutPeakWidth, donutWindow, donutOut);
        if (!passesDonut(donutOut, donutRatioAll, donutRatioOne))
            continue;
        candVal[nRidgeCand] = val;
        candB1[nRidgeCand]  = b1;
        candB2[nRidgeCand]  = b2;
        candPct[nRidgeCand] = pct;
        nRidgeCand++;
    }
}
if (nRidgeCand > 0) {
    candVal = Array.trim(candVal, nRidgeCand);
    order   = Array.rankPositions(candVal);             // ascending
    accB1 = newArray(nRidgeCand);
    accB2 = newArray(nRidgeCand);
    for (oi = nRidgeCand - 1; oi >= 0; oi--) {
        c = order[oi];
        b1 = candB1[c];
        b2 = candB2[c];
        dup = false;
        for (k = 0; k < nRidgeAdded; k++)
            if (abs(b1 - accB1[k]) <= 3 && abs(b2 - accB2[k]) <= 3)
                dup = true;
        if (dup)
            continue;
        accB1[nRidgeAdded] = b1;
        accB2[nRidgeAdded] = b2;
        nRidgeAdded++;
        pct = candPct[c];
        row = nResults;
        setResult("bin1", row, b1);
        setResult("bin2", row, b2);
        setResult("bp1", row, b1 * binSizeBp);
        setResult("bp2", row, b2 * binSizeBp);
        genomeBin1 = undoRecenter(b1, recenterBin, imageSize);
        genomeBin2 = undoRecenter(b2, recenterBin, imageSize);
        setResult("genome_bin1", row, genomeBin1);
        setResult("genome_bin2", row, genomeBin2);
        setResult("genome_bp1", row, genomeBin1 * binSizeBp);
        setResult("genome_bp2", row, genomeBin2 * binSizeBp);
        setResult("n_detections", row, 0);
        setResult("note", row, "ridge region (donut test)");
        setResult("outline_ids", row, "");
        setResult("mean_area", row, NaN);
        setResult("mean_intensity", row, NaN);
        setResult("loop_strength", row, NaN);
        setResult("centroid1", row, NaN);
        setResult("centroid2", row, NaN);
        setResult("refine_shift1", row, 0);
        setResult("refine_shift2", row, 0);
        setResult("log2_obs_exp", row, candVal[c]);
        setResult("oe_percentile", row, pct);
        if (pct >= confidenceHighPercentile)
            setResult("confidence", row, "high");
        else if (pct >= confidenceMediumPercentile)
            setResult("confidence", row, "medium");
        else
            setResult("confidence", row, "low");
    }
}
updateResults();
print("Step 7c ridge search: " + nRidgeMax + " local maxima inside the ridge mask (>= " + minLoopDist +
      " bins from the diagonal); " + nRidgePct + " at or above the " + (100 * ridgeMinPercentile) +
      "th distance-matched percentile; " + nRidgeCand + " passed the donut test; " + nRidgeAdded +
      " added after merging calls within 3 bins.");


// integrated_excess for every loop, main-pass and ridge calls alike.
for (row = 0; row < nResults; row++) {
    b1 = round(getResult("bin1", row));
    b2 = round(getResult("bin2", row));
    sum = 0;
    for (dy = -integratedExcessHalf; dy <= integratedExcessHalf; dy++) {
        for (dx = -integratedExcessHalf; dx <= integratedExcessHalf; dx++) {
            yy = ((b1 + dy) % imageSize + imageSize) % imageSize;
            xx = ((b2 + dx) % imageSize + imageSize) % imageSize;
            sum = sum + excArr[yy * imageSize + xx];
        }
    }
    setResult("integrated_excess", row, 100 * sum / nonLocalMedian);
}
updateResults();
for (row = 0; row < nResults; row++) {             // donut ratios for every loop (a check; only ridge calls must pass)
    donutRatios(round(getResult("bin1", row)), round(getResult("bin2", row)), obsArr, expected, imageSize,
                donutPeakWidth, donutWindow, donutOut);
    setResult("donut_ratio", row, donutOut[0]);
    setResult("lowerleft_ratio", row, donutOut[1]);
    setResult("horizontal_ratio", row, donutOut[2]);
    setResult("vertical_ratio", row, donutOut[3]);
}
updateResults();
print("Step 7c: log2_obs_exp and confidence read from the log2(obs/exp) map with a rolling ball of radius " +
      localBackgroundRadius + " bins subtracted; integrated_excess summed over " +
      (2 * integratedExcessHalf + 1) + " x " + (2 * integratedExcessHalf + 1) +
      " bins, as % of the median bin's non-local contacts (" + nonLocalMedian + ").");

// Display the log2(obs/exp) map on the same red/blue scale as HiC_APA.
selectWindow(obsExpTitle);
applyRedBlueLut();
setMinAndMax(-1, 1);

// ---- 7d. Mark every resolved loop on a copy of Loop_density_smoothed ------
// Each row in the Results table is a resolved loop with an unordered
// (bin1=lo, bin2=hi) anchor pair, but on the actual matrix a loop shows up
// at TWO symmetric pixel positions, (bin1,bin2) and (bin2,bin1). Only one is
// drawn -- column = bin2 (larger), row = bin1 (smaller) -- which is always
// the point above and to the right of the diagonal, consistently, loop
// after loop.
//
// bin1/bin2 are already in Loop_density_smoothed's own (original, not upsampled)
// pixel scale -- see the calibration note in step 7 above -- so they're
// multiplied by loopCircleUpsampleFactor below to land on the right pixel
// in this step's own upsampled working copy.
//
// Each marker is drawn directly into the Overlay of a disposable working
// copy (Overlay.drawEllipse for circles, Overlay.drawLine for the plus signs
// of ridge-search loops), Flattened into the real output image, and the
// disposable copy closed. Drawing into the Overlay instead of adding
// selections to it means no selection -- whose handle squares Flatten would
// bake in -- ever exists. The working copy
// is upsampled first (bilinear) so the circles -- drawn at
// loopCircleUpsampleFactor times the pixel density -- come out smooth
// rather than blocky when the image is zoomed in on.
selectWindow(loopDensitySmoothedTitle);
run("Select None");
run("Duplicate...", "title=" + loopsSourceTitle);
loopCircleImageSize = imageSize * loopCircleUpsampleFactor;
run("Size...", "width=" + loopCircleImageSize + " height=" + loopCircleImageSize + " interpolation=Bilinear");
run("RGB Color");   // bake loopDensitySmoothedTitle's current LUT into real RGB pixels first, so
                     // loopCircleColor renders as an actual color, not through the LUT

// +0.5 before scaling: bin1/bin2 are pixel INDICES (pixel i spans the
// continuous range [i, i+1) in ImageJ's own ROI/overlay coordinate system),
// but a circle's center needs to land on that pixel's visual MIDDLE
// (i+0.5), not its corner (i) -- otherwise every circle is drawn a
// consistent half-bin toward smaller row/column than the peak it's marking.
// This is a different conversion from the floor() used when centroids are
// first turned into bin1/bin2 in resolveSymmetricLoops() -- that one goes
// continuous-coordinate -> bin index; this one goes bin index -> drawing
// coordinate, the opposite direction, so it adds back the 0.5 that floor()
// deliberately drops.
nLoops = nResults;
for (i = 0; i < nLoops; i++) {
    loopBin1 = getResult("bin1", i);
    loopBin2 = getResult("bin2", i);
    cx = (loopBin2 + 0.5) * loopCircleUpsampleFactor;   // column = the larger bin -> right of the diagonal
    cy = (loopBin1 + 0.5) * loopCircleUpsampleFactor;   // row    = the smaller bin -> above the diagonal
    r  = loopCircleRadius * loopCircleUpsampleFactor;
    // Markers are drawn straight into the overlay (Overlay.drawEllipse /
    // Overlay.drawLine) rather than as selections added to it, so no
    // selection -- and none of its handle squares -- ever exists for
    // Flatten to bake into the image.
    setLineWidth(loopCircleStrokeWidth * loopCircleUpsampleFactor);
    if (startsWith(getResultString("note", i), "ridge region")) {
        // ridge-search loop: a plus sign, arms the same length as the circle radius
        setColor(ridgeCallColor);
        Overlay.drawLine(cx - r, cy, cx + r, cy);
        Overlay.drawLine(cx, cy - r, cx, cy + r);
    } else {
        if (getResultString("confidence", i) == "low")
            setColor(lowConfidenceCircleColor);
        else
            setColor(loopCircleColor);
        Overlay.drawEllipse(cx - r, cy - r, r * 2, r * 2);
    }
}
// Outline of the ridge mask (the region the ridge search covers), at this image's scale.
addRidgeOutline(ridgeMaskTitle, ridgeOutlineTitle, loopsSourceTitle, imageSize, loopCircleUpsampleFactor,
                ridgeOutlineColor, ridgeOutlineWidth * loopCircleUpsampleFactor);
selectWindow(loopsSourceTitle);
Overlay.show;
run("Select None");     // belt-and-suspenders: nothing selected when Flatten runs
run("Flatten");
rename(loopsTitle);
close(loopsSourceTitle);
// ... and on Loop_density_smoothed itself, as an overlay (kept when saved as TIFF)
addRidgeOutline(ridgeMaskTitle, ridgeOutlineTitle, loopDensitySmoothedTitle, imageSize, 1,
                ridgeOutlineColor, ridgeOutlineWidth);
close(ridgeMaskTitle);

// ---- 7e. Aggregate Peak Analysis (APA): log2(observed/expected) around every loop --
// Averages a small square of the Hi-C map, centered on every resolved loop,
// into one image -- the standard "aggregate peak analysis" check: real loops
// produce a focused enrichment in the middle of the average, while
// noise/false positives wash out since they aren't consistently centered.
//
// Built from mainTitle (HiC_map: the scaled, diagonal-smoothed contact map,
// before any background subtraction or ridge removal) rather than from
// Loop_density, so the result is directly comparable to APA/pileup plots
// from other loop callers. Expressed the same way coolpup.py/cooltools
// pileups are: every pixel is divided by the EXPECTED contact frequency at
// its genomic distance, the ratios are averaged across loops, and the
// average is shown as log2 -- so 0 means "no more contact than expected at
// that distance", positive (red) is enrichment and negative (blue) is
// depletion.
//
// Expected: the mean of every HiC_map pixel at the same CIRCULAR genomic
// distance, computed once in step 7c. Pixels <= 0 are treated as missing
// (bins ICE filtered out) and skipped in the average below.
//
// mainTitle is still at step 2's PADDED canvas size, wrapped periodically, so
// a loop's position in it is (bin + padding/2), and a window near the image
// edge reads genuinely wrapped signal rather than empty padding. bin1/bin2
// are in the UNPADDED pixel scale (0..imageSize-1) and can land on a
// half-integer (see resolveSymmetricLoops()'s rounding-boundary case), so
// each is rounded to the nearest whole pixel before use. Uses column = bin2
// (larger), row = bin1 (smaller) -- the same single upper-right-of-diagonal
// position step 7d draws its loop circles at.
// Uses the same expected[] array step 7c built.
apaHalf  = (apaSize - 1) / 2;

selectWindow(mainTitle);
run("Select None");

// Average obs/expected per window pixel, across loops. apaPixCount is kept
// per pixel since missing (<= 0) pixels are skipped individually.
apaSum      = newArray(apaSize * apaSize);   // row-major accumulators, zero-filled by default
apaPixCount = newArray(apaSize * apaSize);
nLoopsForApa = nResults;
apaCount = 0;
for (i = 0; i < nLoopsForApa; i++) {
    loopRow = round(getResult("bin1", i));
    loopCol = round(getResult("bin2", i));
    for (dy = -apaHalf; dy <= apaHalf; dy++) {
        for (dx = -apaHalf; dx <= apaHalf; dx++) {
            v = getPixel(loopCol + dx + padHalf, loopRow + dy + padHalf);
            d = abs((loopCol + dx) - (loopRow + dy)) % imageSize;
            d = minOf(d, imageSize - d);
            if (v > 0 && expected[d] > 0) {
                idx = (dy + apaHalf) * apaSize + (dx + apaHalf);
                apaSum[idx]      = apaSum[idx] + v / expected[d];
                apaPixCount[idx] = apaPixCount[idx] + 1;
            }
        }
    }
    apaCount++;
}

newImage(apaTitle, "32-bit black", apaSize, apaSize, 1);
selectWindow(apaTitle);
for (dy = 0; dy < apaSize; dy++) {
    for (dx = 0; dx < apaSize; dx++) {
        idx = dy * apaSize + dx;
        if (apaPixCount[idx] > 0)
            setPixel(dx, dy, log(apaSum[idx] / apaPixCount[idx]) / log(2));
    }
}
apaCenterLog2 = getPixel(apaHalf, apaHalf);

// Red/blue LUT, white at 0, with a display range symmetric around 0 --
// either fixed (apaColorLimit > 0, so different conditions share one color
// scale) or, if apaColorLimit is 0, the largest |value| in this APA.
applyRedBlueLut();
if (apaColorLimit > 0) {
    apaLimit = apaColorLimit;
} else {
    getStatistics(area, mean, apaMin, apaMax);
    apaLimit = maxOf(abs(apaMin), abs(apaMax));
    if (apaLimit == 0)
        apaLimit = 1;   // flat/empty APA -- avoid a zero-width display range
}

// Enlarge for viewing. interpolation=None keeps every averaged value exactly
// as computed -- each one just becomes an apaUpsampleFactor x
// apaUpsampleFactor block -- so this window still holds the real 32-bit
// log2(obs/exp) values, only easier to see.
apaDisplaySize = apaSize * apaUpsampleFactor;
run("Size...", "width=" + apaDisplaySize + " height=" + apaDisplaySize + " depth=1 interpolation=None");
setMinAndMax(-apaLimit, apaLimit);
print("Step 7e APA: averaged " + apaCount + " loop(s) into a " + apaSize + "x" + apaSize +
      " log2(obs/exp) image (" + apaTitle + "), enlarged to " + apaDisplaySize + "x" +
      apaDisplaySize + " px. Center pixel log2(obs/exp) = " + apaCenterLog2 +
      "; color range +/-" + apaLimit + ".");

// ---- 8. Integrated_loop_density, for comparison to ChIP-seq / RNA-seq ----
selectWindow(obsExpTitle);   // cropped here too, now that step 7e no longer needs the padding
run("Canvas Size...", "width=" + imageSize + " height=" + imageSize + " position=Center");

selectWindow(loopDensityTitle);
run("Canvas Size...", "width=" + imageSize + " height=" + imageSize + " position=Center");
setMinAndMax(0, displayMax);
makeRectangle(0, 0, imageSize, imageSize);   // explicit rectangle rather than "Select All"

// Built via getProfile() -- the same column-average track run("Plot Profile")
// would show -- rather than letting Plot Profile build its own plot, so the
// x-axis can be real genome position in bp instead of raw TIFF pixel index.
profile = getProfile();
nProfile = profile.length;
xGenomeBp = newArray(nProfile);
for (i = 0; i < nProfile; i++)
    xGenomeBp[i] = undoRecenterSigned(i, recenterBin, imageSize) * binSizeBp;

Plot.create(integratedLoopDensityTitle, "Genome position (bp)", "ILD (a.u.)", xGenomeBp, profile);
Plot.show();

// ---- 9. Save outputs -----------------------------------------------------
// Saves six deliverables to saveDir, each named with imageIdentifier
// pasted onto the end -- skipped entirely (with a note to the Log) if no
// save directory was given in the step 0 dialog. All windows/tables stay
// open afterward too, for inspection or a manual re-save.
if (saveDir == "") {
    print("No save directory given in the startup dialog -- skipping automatic file saving.");
} else {
    // saveAs("Tiff", ...) renames the active window to match the saved
    // filename (e.g. loopsTitle becomes "Loop_density_smoothed_loops_<id>")
    // as a side effect -- so without the rename() right after each one below,
    // the window would keep that renamed title for the rest of the ImageJ
    // session, and the closeIfOpen(loopsTitle)/etc. calls at the very
    // start of the macro (which only know the canonical, un-suffixed titles)
    // would no longer find it on a later run, leaking a stale window instead
    // of clearing it. Renaming back immediately keeps every window under its
    // canonical title regardless of whether saving happened on a previous
    // run, so startup cleanup always works the same way.
    // Loop_density is saved 16-bit, cropped back to the original image size
    // (step 8), with the step 8 selection cleared so it isn't stored in the file.
    selectWindow(loopDensityTitle);
    run("Select None");
    saveAs("Tiff", saveDir + loopDensityTitle + "_" + imageIdentifier + ".tif");
    rename(loopDensityTitle);

    selectWindow(loopsTitle);
    saveAs("Tiff", saveDir + loopsTitle + "_" + imageIdentifier + ".tif");
    rename(loopsTitle);

    selectWindow(loopDensitySmoothedTitle);
    saveAs("Tiff", saveDir + loopDensitySmoothedTitle + "_" + imageIdentifier + ".tif");
    rename(loopDensitySmoothedTitle);

    // HiC_APA is saved as a PNG for easy viewing: a disposable copy is
    // converted to RGB first, which bakes in the red/blue LUT and the current
    // display range (ScaleConversions is back on by this point -- see step 5),
    // and closed once saved. apaTitle itself stays open, still 32-bit, for
    // reading off the actual averaged values.
    selectWindow(apaTitle);
    run("Select None");
    run("Duplicate...", "title=" + apaTitle + "_png");
    run("RGB Color");
    saveAs("PNG", saveDir + apaTitle + "_" + imageIdentifier + ".png");
    close();

    saveAs("Results", saveDir + "loop_calls_" + imageIdentifier + ".csv");

    // The Integrated_loop_density plot data, written directly from the xGenomeBp/profile
    // arrays already built above -- simpler and more reliable than reading the
    // values back out of the Plot window.
    plotCsvPath = saveDir + "integrated_loop_density_" + imageIdentifier + ".csv";
    plotCsv = File.open(plotCsvPath);
    print(plotCsv, "genome_position_bp,ILD_au");
    for (i = 0; i < xGenomeBp.length; i++)
        print(plotCsv, xGenomeBp[i] + "," + profile[i]);
    File.close(plotCsv);

    print("Saved outputs to " + saveDir + " with identifier \"" + imageIdentifier + "\":");
    print("  " + loopDensityTitle + "_" + imageIdentifier + ".tif");
    print("  " + loopsTitle + "_" + imageIdentifier + ".tif");
    print("  " + loopDensitySmoothedTitle + "_" + imageIdentifier + ".tif");
    print("  " + apaTitle + "_" + imageIdentifier + ".png");
    print("  loop_calls_" + imageIdentifier + ".csv");
    print("  integrated_loop_density_" + imageIdentifier + ".csv");
}
setBatchMode("exit and display");   // show every window left open, now that processing is done

// ============================================================
//  Functions
// ============================================================

// Replaces every pixel on the main diagonal (x == y) of the CURRENT image
// with the median of its 3x3 neighborhood. Works by running the median
// filter on a disposable full-image duplicate, then copying just the
// diagonal pixel values back onto the real image.
function smoothDiagonal(title, size) {
    // The 3x3 neighbourhood wraps around the circular map, so the two corner pixels (0, 0) and
    // (size-1, size-1) are treated like every other diagonal pixel and the result does not depend
    // on where the map starts.
    selectWindow(title);
    run("Select None");
    diagValues = newArray(size);
    nb = newArray(9);
    for (x = 0; x < size; x++) {
        k = 0;
        for (dy = -1; dy <= 1; dy++)
            for (dx = -1; dx <= 1; dx++) {
                nb[k] = getPixel((x + dx + size) % size, (x + dy + size) % size);
                k++;
            }
        Array.sort(nb);
        diagValues[x] = nb[4];
    }
    for (x = 0; x < size; x++)
        setPixel(x, x, diagValues[x]);
}

// Returns a robust estimate of the CURRENT image's ICE row/column-sum
// constant: the median of every row's pixel sum. The median is used rather
// than the mean because a minority of rows/columns correspond to bins ICE
// filtered out before balancing -- those rows sum to ~0 and would drag a
// mean down. Row sums are read via getStatistics on a 1-row selection
// rather than a getPixel loop, since the per-pixel summation then runs in
// ImageJ's compiled statistics code instead of the macro interpreter.
function computeRowSum(title) {
    selectWindow(title);
    run("Select None");
    width  = getWidth();
    height = getHeight();
    rowSums = newArray(height);
    for (y = 0; y < height; y++) {
        makeRectangle(0, y, width, 1);
        getStatistics(area, mean, min, max, std, histogram);
        rowSums[y] = mean * width;
    }
    run("Select None");
    ranks = Array.rankPositions(rowSums);
    return rowSums[ranks[floor(height / 2)]];   // middle-ranked row sum -- an adequate median for this
}

// Returns the multiplier that puts the `percentile`-th percentile (0-1) of
// the CURRENT image's raw pixel values at gray level `targetLevel`. A high
// percentile (e.g. 0.999) rather than the true max is used as the scaling
// ceiling because Hi-C matrices routinely have a handful of extreme outlier
// pixels that would otherwise single-handedly determine the scale.
//
// getHistogram's 5-argument form is required for 16/32-bit images, and its
// "values" output isn't reliable on every ImageJ installation, so the
// percentile's actual value is reconstructed from the bin index and the
// known bin width instead.
function computeScaleFactor(title, percentile, targetLevel) {
    selectWindow(title);
    getStatistics(area, mean, minVal, maxVal, std, histogram);
    if (maxVal <= 0)
        return 1;   // degenerate/empty image -- avoid a divide-by-zero below

    nBins = 10000;   // fine enough to locate the percentile precisely without needing "values"
    getHistogram(values, counts, nBins, minVal, maxVal);

    total = 0;
    for (i = 0; i < counts.length; i++)
        total += counts[i];

    target = percentile * total;
    cum = 0;
    binIndex = counts.length - 1;   // fallback: top bin, in case rounding leaves the loop short
    for (i = 0; i < counts.length; i++) {
        cum += counts[i];
        if (cum >= target) {
            binIndex = i;
            break;
        }
    }

    binWidth = (maxVal - minVal) / nBins;
    percentileValue = minVal + binIndex * binWidth;
    if (percentileValue <= 0)
        percentileValue = maxVal;   // guard against landing in an all-zero bin

    return targetLevel / percentileValue;
}

// Computes robust background statistics for the CURRENT image -- the median
// and an outlier-resistant stand-in for the standard deviation (robustSD =
// IQR / 1.349, from the interquartile range) -- and returns them to the
// caller as a single string "robustMedian,robustSD" (split() + parseFloat()
// on the caller's side turns it back into two numbers).
//
// IMPORTANT ImageJ macro quirk this works around: every variable a function
// assigns internally is local to that function, even when a same-named
// variable already exists at the macro's top level -- pre-declaring
// robustMedian/robustSD as top-level globals before calling this function
// would NOT make an assignment inside the function visible outside it. A
// single return value is the only reliable mechanism here, so both numbers
// are packed into one comma-delimited string and returned that way instead.
// See the "Known ImageJ macro quirks" section of the README for more on this.
//
// A second, related quirk applies to building that return string: writing it
// as `robustMedian + "," + robustSD` (a NUMBER first) evaluates to the
// numeric value NaN instead of a string, because this macro dialect's "+"
// does plain arithmetic whenever the LEFT operand of a given "+" is a
// number -- it does not fall back to concatenation just because the right
// operand is a string, the way STRING + NUMBER does everywhere else in this
// file (e.g. every `print("label=" + someNumber)` call, which always leads
// with a string). A leading `""` forces concatenation mode from the start:
// `"" + robustMedian + "," + robustSD`. Also in the "Known ImageJ macro
// quirks" section of the README.
//
// Works for any bit depth, INCLUDING a 32-bit image with negative values --
// unlike a fixed 0-65535 range, it bins over the image's own actual
// [min, max] (via getStatistics), the same approach computeScaleFactor()
// above uses. Q1, median, and Q3 are found in one pass over the histogram,
// in one function, since macro arrays don't reliably survive being passed
// between user-defined functions.
function computeRobustStats(title) {
    selectWindow(title);
    getStatistics(area, mean, minVal, maxVal, std, histogram);
    if (maxVal <= minVal) {
        // degenerate/constant image -- avoid a divide-by-zero below
        return "" + minVal + ",0";
    }

    nBins = 65536;   // fine enough to locate the quartiles precisely without needing "values"
    getHistogram(values, counts, nBins, minVal, maxVal);
    binWidth = (maxVal - minVal) / nBins;

    total = 0;
    for (i = 0; i < counts.length; i++)
        total += counts[i];

    q1Target     = 0.25 * total;
    medianTarget = 0.50 * total;
    q3Target     = 0.75 * total;
    q1Bin = -1; medianBin = -1; q3Bin = -1;
    cum = 0;
    for (i = 0; i < counts.length; i++) {
        cum += counts[i];
        if (q1Bin < 0 && cum >= q1Target)
            q1Bin = i;
        if (medianBin < 0 && cum >= medianTarget)
            medianBin = i;
        if (q3Bin < 0 && cum >= q3Target) {
            q3Bin = i;
            break;
        }
    }

    q1           = minVal + q1Bin * binWidth;
    robustMedian = minVal + medianBin * binWidth;
    q3           = minVal + q3Bin * binWidth;
    robustSD     = (q3 - q1) / 1.349;
    return "" + robustMedian + "," + robustSD;
}

// Sets a threshold on the CURRENT image at baselineValue + sdMultiplier * sd,
// selecting pixels AT OR ABOVE that value as foreground. baselineValue and sd
// are supplied by the caller rather than computed here, so this is just a
// thin wrapper around setThreshold() that does the arithmetic. See the step 7
// call site for why its sd is one measured earlier, elsewhere, rather than
// fresh from the current image's own (possibly skewed) histogram.
function setRobustThreshold(sdMultiplier, baselineValue, sd) {
    getStatistics(area, mean, minVal, maxVal, std, histogram);
    thresholdValue = baselineValue + sdMultiplier * sd;
    setThreshold(thresholdValue, maxVal);
}

// Maps bin1/bin2 -- a 0-based pixel/bin index into THIS TIFF -- back to the
// REAL genome bin it came from before hic_circular_recenter_to_tiff.py's
// --recenter-bin rotation was applied. Returns tiffBin unchanged if
// recenterBin is 0 (no recentering was used).
//
// hic_circular_recenter_to_tiff.py builds the recentered matrix, in 0-based
// numpy terms, as recentered[r, c] = original[(start + r) mod genomeBins, ...],
// start = recenterBin - 1 (recenterBin is the script's 1-based CLI argument).
// Since tiffBin IS that same 0-based r, the original 0-based index is just
// (start + tiffBin) mod genomeBins. recenterBin >= 1 and tiffBin >= 0, so
// recenterBin - 1 + tiffBin is always >= 0 -- no negative-operand modulo edge case.
function undoRecenter(tiffBin, recenterBin, genomeBins) {
    if (recenterBin == 0)
        return tiffBin;
    return (recenterBin - 1 + tiffBin) % genomeBins;
}

// Same tiffBin -> real genome bin mapping as undoRecenter(), but for
// building a CONTINUOUS axis across a whole recentered TIFF (step 8's
// Integrated_loop_density plot) rather than looking up one discrete loop
// call at a time. undoRecenter() always returns a value in [0, genomeBins),
// which is correct for a single bin, but calling it for every tiffBin in
// sequence and plotting the results produces a sawtooth: the track counts
// up toward genomeBins, then drops straight back down to 0 right where the
// recentering wrapped the genome around, even though nothing discontinuous
// actually happens there biologically.
//
// This instead lets the segment that appears BEFORE that wrap point in the
// TIFF come out NEGATIVE -- that many bins before real genome bin 0,
// wherever genome bin 0 happens to land in this TIFF -- so the whole track
// increases monotonically with no jump in the middle.
//
// If recenterBin is 0 (no recentering) or 1 (recentering was a no-op) there's
// no wrap to avoid, so this returns the same non-negative value undoRecenter() would.
function undoRecenterSigned(tiffBin, recenterBin, genomeBins) {
    if (recenterBin == 0)
        return tiffBin;
    rawBin = (recenterBin - 1) + tiffBin;
    genomeBin = rawBin % genomeBins;
    if (rawBin < genomeBins && recenterBin > 1)
        genomeBin = genomeBin - genomeBins;
    return genomeBin;
}

// Collapses the Results table's raw particle detections (one row per
// particle found on EITHER side of the diagonal) into one row per loop, and
// rewrites the Results table in place with that collapsed version.
//
// Every row's (X,Y) centroid is converted to a bin and reduced to an
// unordered (lo,hi) pair (lo = the smaller of the two, hi = the larger),
// since a loop at column A / row B and the same loop at column B / row A
// are the same loop. From there, three cases can happen:
//
//   1. Two (or more) rows land on the exact same (lo,hi) bin pair -- the
//      expected mirror case, where the same loop was found once above and
//      once below the diagonal at the same bin. Averaged into one row.
//
//   2. Two rows differ by exactly 1 bin in ONE of the two coordinates, but
//      not both -- a loop whose true centroid sits right on a rounding
//      boundary, so the two sides' scans round it in opposite directions.
//      Found with a second, wider pass over whatever case 1 left unmatched,
//      and averaged into one row the same way.
//
//   3. A row that doesn't pair up under either case above -- a genuine
//      single-sided detection. Kept as its own row, flagged via the "note"
//      column, so it's easy to go back and check whether it's a real
//      weak/borderline call or a fluke.
//
// Both matching cases (1 and 2) additionally require the two rows to come
// from OPPOSITE sides of the diagonal -- see the `side` array below. Without
// that check, two distinct same-side detections that happen to sit close
// together (e.g. a watershed split, or two genuinely separate weak calls)
// could be mistaken for a mirror pair of each other, which would silently
// merge them into one row and strand the real mirror detection as a
// spurious single-sided one instead.
function resolveSymmetricLoops(genomeBins, recenterBin, binSizeBp) {
    n = nResults;
    if (n == 0) {
        print("resolveSymmetricLoops: Results table is empty, nothing to do.");
        return;
    }

    // ---- Read the raw particle table before it gets rewritten ----
    rawX    = newArray(n);
    rawY    = newArray(n);
    rawArea = newArray(n);
    rawMean = newArray(n);
    for (i = 0; i < n; i++) {
        rawX[i]    = getResult("X", i);
        rawY[i]    = getResult("Y", i);
        rawArea[i] = getResult("Area", i);
        rawMean[i] = getResult("Mean", i);
    }
    // Unrounded centroids, lower/upper coordinate, kept for the centroid1/centroid2
    // columns (averaged over a mirror pair the same way the bins are).
    rawLo = newArray(n);
    rawHi = newArray(n);
    for (i = 0; i < n; i++) {
        rawLo[i] = minOf(rawX[i], rawY[i]);
        rawHi[i] = maxOf(rawX[i], rawY[i]);
    }

    // ---- Convert centroids to bins, and build a single sortable key per row ----
    // ImageJ's particle centroids are in the convention where pixel index i
    // occupies the continuous range [i, i+1) -- so a particle sitting
    // squarely on pixel 245 is reported as X/Y ~= 245.5, not 245. floor() is
    // the correct conversion back to a pixel/bin index -- round() would bump
    // the bin up by 1 whenever the fractional part is >= 0.5, even though
    // the centroid is still physically inside the lower-numbered pixel.
    //
    // key = lo*100000 + hi. This assumes hi never reaches 100000 bins,
    // which holds for any Hi-C map this macro would realistically be run
    // on (100000 bins at 10kb would be a 1 Gb genome).
    //
    // side records which side of the diagonal each raw detection came from:
    // +1 for row < column (Y < X) -- the "above and to the right" position
    // step 7d draws its loop circles at -- -1 for the mirror position below
    // the diagonal (Y > X), and 0 if the detection's bin sits exactly ON the
    // diagonal (bx == by; ambiguous, so treated as compatible with either
    // side in the matching passes below rather than excluded from both).
    binLo = newArray(n);
    binHi = newArray(n);
    side  = newArray(n);
    key   = newArray(n);
    for (i = 0; i < n; i++) {
        bx = floor(rawX[i]);
        by = floor(rawY[i]);
        binLo[i] = minOf(bx, by);
        binHi[i] = maxOf(bx, by);
        key[i]   = binLo[i] * 100000 + binHi[i];
        if (by < bx)
            side[i] = 1;
        else if (by > bx)
            side[i] = -1;
        else
            side[i] = 0;
    }

    // Array.rankPositions(key) returns the indexes that sort key ascending
    // -- i.e. order[0] is the index of the smallest key, order[1] the next
    // smallest, and so on -- which is what lets identical/adjacent keys end
    // up next to each other in the pass below.
    order = Array.rankPositions(key);

    // Output rows, pre-allocated at the input size -- an upper bound on how
    // many rows the collapsed table could possibly need, since collapsing
    // can only ever reduce the row count. outCount tracks how many of these
    // pre-allocated slots are actually filled, since macro arrays can't
    // grow on their own.
    outBin1  = newArray(n);
    outBin2  = newArray(n);
    outArea  = newArray(n);
    outMean  = newArray(n);
    outNDet  = newArray(n);
    outNote  = newArray(n);   // 0 = mirror (clean or rounding-boundary), 1 = true singleton
    outIds   = newArray(n);   // comma-separated outline_id(s) this row came from
    outC1    = newArray(n);   // unrounded centroid, lower coordinate (TIFF bins)
    outC2    = newArray(n);   // unrounded centroid, upper coordinate
    outCount = 0;

    // Rows that don't share an exact key with anything else -- or that do,
    // but not with anything from the opposite side (see below) -- set aside
    // for the wider rounding-boundary search in pass 2. pendId records each
    // pending row's outline_id: order[k] is the original raw-Results row
    // index, which (per the "add" + sequential rename in step 7) is exactly
    // that particle's ROI Manager slot / burned-in label number.
    pendLo    = newArray(n);
    pendHi    = newArray(n);
    pendArea  = newArray(n);
    pendMean  = newArray(n);
    pendSide  = newArray(n);
    pendId    = newArray(n);
    pendC1    = newArray(n);
    pendC2    = newArray(n);
    pendCount = 0;

    // ---- Pass 1: collapse exact-key groups (case 1 above) ----
    // Deliberately written with nested if/break rather than a single
    // compound "j + 1 < n && key[order[j+1]] == ..." while-condition --
    // && does not short-circuit on every ImageJ installation, so a combined
    // condition like that can still evaluate key[order[j+1]] out of bounds.
    i = 0;
    while (i < n) {
        j = i;
        while (true) {
            if (j + 1 >= n)
                break;
            if (key[order[j + 1]] != key[order[i]])
                break;
            j++;
        }
        groupSize = j - i + 1;

        // A same-key group only counts as a genuine mirror match if it has
        // detections from both sides of the diagonal (an on-diagonal, side-0
        // detection is compatible with either) -- otherwise every member is
        // on the same side, which isn't a real mirror pair, so each is
        // pushed to pending individually instead of being averaged together.
        sawPlus = false; sawMinus = false; sawZero = false;
        for (k = i; k <= j; k++) {
            if (side[order[k]] == 1)
                sawPlus = true;
            else if (side[order[k]] == -1)
                sawMinus = true;
            else
                sawZero = true;
        }
        crossSide = sawZero || (sawPlus && sawMinus);

        if (groupSize >= 2 && crossSide) {
            sumArea = 0;
            sumMean = 0;
            sumC1   = 0;
            sumC2   = 0;
            idStr   = "" + order[i];
            for (k = i; k <= j; k++) {
                sumArea += rawArea[order[k]];
                sumMean += rawMean[order[k]];
                sumC1   += rawLo[order[k]];
                sumC2   += rawHi[order[k]];
                if (k > i)
                    idStr = idStr + "," + order[k];
            }
            outBin1[outCount] = binLo[order[i]];
            outBin2[outCount] = binHi[order[i]];
            outArea[outCount] = sumArea / groupSize;
            outMean[outCount] = sumMean / groupSize;
            outNDet[outCount] = groupSize;
            outNote[outCount] = 0;
            outIds[outCount]  = idStr;
            outC1[outCount]   = sumC1 / groupSize;
            outC2[outCount]   = sumC2 / groupSize;
            outCount++;
        } else {
            for (k = i; k <= j; k++) {
                pendLo[pendCount]   = binLo[order[k]];
                pendHi[pendCount]   = binHi[order[k]];
                pendArea[pendCount] = rawArea[order[k]];
                pendMean[pendCount] = rawMean[order[k]];
                pendSide[pendCount] = side[order[k]];
                pendId[pendCount]   = order[k];
                pendC1[pendCount]   = rawLo[order[k]];
                pendC2[pendCount]   = rawHi[order[k]];
                pendCount++;
            }
        }
        i = j + 1;
    }

    // ---- Pass 2: pair up rounding-boundary singletons (case 2 above) ----
    // pendCount is small in practice, so the plain O(pendCount^2) search
    // here is in no way a performance concern.
    used = newArray(pendCount);   // 0/1 flags; newArray() zero-fills by default
    for (a = 0; a < pendCount; a++) {
        if (used[a] == 1)
            continue;
        for (b = a + 1; b < pendCount; b++) {
            if (used[b] == 1)
                continue;
            // Same-side pair (both nonzero and equal) -- not a real mirror
            // candidate, however close their bins are; keep searching rather
            // than letting two same-side detections match each other.
            if (pendSide[a] != 0 && pendSide[b] != 0 && pendSide[a] == pendSide[b])
                continue;
            sameLo = (pendLo[a] == pendLo[b]) && (abs(pendHi[a] - pendHi[b]) == 1);
            sameHi = (pendHi[a] == pendHi[b]) && (abs(pendLo[a] - pendLo[b]) == 1);
            if (sameLo || sameHi) {
                outBin1[outCount] = (pendLo[a] + pendLo[b]) / 2;
                outBin2[outCount] = (pendHi[a] + pendHi[b]) / 2;
                outArea[outCount] = (pendArea[a] + pendArea[b]) / 2;
                outMean[outCount] = (pendMean[a] + pendMean[b]) / 2;
                outNDet[outCount] = 2;
                outNote[outCount] = 0;
                outIds[outCount]  = "" + pendId[a] + "," + pendId[b];
                outC1[outCount]   = (pendC1[a] + pendC1[b]) / 2;
                outC2[outCount]   = (pendC2[a] + pendC2[b]) / 2;
                outCount++;
                used[a] = 1;
                used[b] = 1;
                b = pendCount;   // break out of the inner loop
            }
        }
    }

    // ---- Pass 3: whatever's left is a true single-sided detection ----
    nSingletons = 0;
    for (a = 0; a < pendCount; a++) {
        if (used[a] == 0) {
            outBin1[outCount] = pendLo[a];
            outBin2[outCount] = pendHi[a];
            outArea[outCount] = pendArea[a];
            outMean[outCount] = pendMean[a];
            outNDet[outCount] = 1;
            outNote[outCount] = 1;
            outIds[outCount]  = "" + pendId[a];
            outC1[outCount]   = pendC1[a];
            outC2[outCount]   = pendC2[a];
            outCount++;
            nSingletons++;
        }
    }

    print("resolveSymmetricLoops: " + n + " raw detections -> " + outCount +
          " unique loops (" + nSingletons + " single-sided).");

    // ---- Rewrite the Results table in place with the collapsed version ----
    // "Clear Results" (not a plain row-delete) is needed to fully reset the
    // table's underlying measurement state -- otherwise every rewritten row
    // can end up showing stale Area/Mean/X/Y values from the last raw particle.
    run("Clear Results");

    for (r = 0; r < outCount; r++) {
        setResult("bin1", r, outBin1[r]);
        setResult("bin2", r, outBin2[r]);
        setResult("bp1", r, outBin1[r] * binSizeBp);
        setResult("bp2", r, outBin2[r] * binSizeBp);
        genomeBin1 = undoRecenter(outBin1[r], recenterBin, genomeBins);
        genomeBin2 = undoRecenter(outBin2[r], recenterBin, genomeBins);
        setResult("genome_bin1", r, genomeBin1);
        setResult("genome_bin2", r, genomeBin2);
        setResult("genome_bp1", r, genomeBin1 * binSizeBp);
        setResult("genome_bp2", r, genomeBin2 * binSizeBp);
        setResult("n_detections", r, outNDet[r]);
        if (outNote[r] == 1)
            setResult("note", r, "TRUE single-sided detection");
        else
            setResult("note", r, "mirror (clean or rounding-boundary)");
        setResult("outline_ids", r, outIds[r]);
        setResult("mean_area", r, outArea[r]);
        setResult("mean_intensity", r, outMean[r]);
        setResult("loop_strength", r, outArea[r] * outMean[r]);
        setResult("centroid1", r, outC1[r]);
        setResult("centroid2", r, outC2[r]);
    }
    updateResults();
}

// Reads an image's unpadded imageSize x imageSize core (bin b at padded pixel
// b + padHalf) into a row-major array, selecting the image by title.
function readUnpadded(title, N, padHalf) {
    selectWindow(title);
    arr = newArray(N * N);
    for (r = 0; r < N; r++)
        for (c = 0; c < N; c++)
            arr[r * N + c] = getPixel(c + padHalf, r + padHalf);
    return arr;
}

// HiCCUPS-style local backgrounds for pixel (i, j) (either orientation): fills out[0..3] with observed /
// locally expected for the donut, lower-left, horizontal and vertical kernels. Each
// kernel's expected = (observed summed over the kernel / distance-expected summed over
// the kernel) x distance-expected at (i, j). Pixels with no data and pixels < 2 bins
// from the diagonal are left out; NaN if the centre has no data or a kernel is empty.
function donutRatios(i, j, obsArr, expected, N, p, w, out) {
    // Orient the pixel along the short arc of the circular chromosome (j lies at most N/2 bins
    // downstream of i), so the lower-left kernel always sits between the loop and the diagonal
    // and the horizontal/vertical kernels do not swap when the map is rotated. A pixel exactly
    // half the chromosome away has two equally short arcs; it gets the mean of both orientations.
    fwd = ((j - i) % N + N) % N;
    if (2 * fwd == N) {
        other = newArray(4);
        donutRatiosOriented(i, j, obsArr, expected, N, p, w, out);
        donutRatiosOriented(j, i, obsArr, expected, N, p, w, other);
        for (k = 0; k < 4; k++)
            out[k] = (out[k] + other[k]) / 2;
    } else if (fwd > N - fwd)
        donutRatiosOriented(j, i, obsArr, expected, N, p, w, out);
    else
        donutRatiosOriented(i, j, obsArr, expected, N, p, w, out);
}

function donutRatiosOriented(i, j, obsArr, expected, N, p, w, out) {
    ci = (i % N + N) % N; cj = (j % N + N) % N;
    centre = obsArr[ci * N + cj];
    dc = abs(ci - cj); dc = minOf(dc, N - dc);
    sumO = newArray(4); sumE = newArray(4);
    for (di = -w; di <= w; di++) {
        for (dj = -w; dj <= w; dj++) {
            r = ((i + di) % N + N) % N; c = ((j + dj) % N + N) % N;
            o = obsArr[r * N + c];
            d = abs(r - c); d = minOf(d, N - d);
            if (!(o > 0) || d < 2)
                continue;
            peak = abs(di) <= p && abs(dj) <= p;
            e = expected[d];
            if (!peak && di != 0 && dj != 0)   { sumO[0] += o; sumE[0] += e; }
            if (di > 0 && dj < 0 && !peak)     { sumO[1] += o; sumE[1] += e; }
            if (abs(di) <= 1 && abs(dj) > p)   { sumO[2] += o; sumE[2] += e; }
            if (abs(dj) <= 1 && abs(di) > p)   { sumO[3] += o; sumE[3] += e; }
        }
    }
    for (k = 0; k < 4; k++) {
        if (!(centre > 0) || sumE[k] <= 0)
            out[k] = NaN;
        else
            out[k] = centre / (sumO[k] / sumE[k] * expected[dc]);
    }
}

// The donut test on donutRatios() output: every ratio >= all, and donut or lower-left >= one.
function passesDonut(r, all, one) {
    for (k = 0; k < 4; k++)
        if (isNaN(r[k]) || r[k] < all)
            return false;
    if (maxOf(r[0], r[1]) >= one)     // not "return maxOf(...) >= one": the macro language would return
        return true;                  // maxOf()'s value and silently drop the comparison
    return false;
}

// Outlines the ridge mask (0 inside) on targetTitle, scaled up by `factor`, as an overlay.
function addRidgeOutline(maskTitle, workTitle, targetTitle, N, factor, color, width) {
    selectWindow(maskTitle);
    run("Select None");
    run("Duplicate...", "title=" + workTitle);
    run("Canvas Size...", "width=" + N + " height=" + N + " position=Center");
    if (factor != 1)
        run("Size...", "width=" + (N * factor) + " height=" + (N * factor) + " interpolation=None");
    setThreshold(0, 0.5);
    run("Create Selection");
    selectWindow(targetTitle);
    run("Restore Selection");
    Overlay.addSelection(color, width);
    run("Select None");
    close(workTitle);
}

// Position refinement pick (step 7c): the offset (dy, dx), within `radius` bins
// of (b1, b2), of the highest value in arr (a imageSize x imageSize map stored
// row-major), skipping offsets that would bring the pair closer than minDist to
// the diagonal. Returned packed as (dy + 500) * 1000 + (dx + 500), since macro
// functions return a single value; ties keep the unmoved pixel.
function refinePick(b1, b2, arr, N, radius, minDist) {
    best = arr[((b1 % N + N) % N) * N + ((b2 % N + N) % N)];
    s1 = 0; s2 = 0;
    for (dy = -radius; dy <= radius; dy++) {
        for (dx = -radius; dx <= radius; dx++) {
            if (circDist(b1 + dy, b2 + dx, N) < minDist)
                continue;
            v = arr[(((b1 + dy) % N + N) % N) * N + (((b2 + dx) % N + N) % N)];
            if (v > best) {
                best = v; s1 = dy; s2 = dx;
            }
        }
    }
    return (s1 + 500) * 1000 + (s2 + 500);
}

// Applies a diverging blue-white-red LUT to the CURRENT image: blue at the
// bottom of the display range, white in the middle, red at the top. Entries
// 127 and 128 are both white, so with a display range symmetric around 0,
// a value of 0 renders white. End colors are ColorBrewer RdBu's.
function applyRedBlueLut() {
    reds   = newArray(256);
    greens = newArray(256);
    blues  = newArray(256);
    for (i = 0; i < 128; i++) {
        t = i / 127;                  // 0 = blue end, 1 = white
        reds[i]   = round(33  + t * (255 - 33));
        greens[i] = round(102 + t * (255 - 102));
        blues[i]  = round(172 + t * (255 - 172));
    }
    for (i = 128; i < 256; i++) {
        t = (i - 128) / 127;          // 0 = white, 1 = red end
        reds[i]   = round(255 + t * (178 - 255));
        greens[i] = round(255 + t * (24  - 255));
        blues[i]  = round(255 + t * (43  - 255));
    }
    setLut(reds, greens, blues);
}

// Circular distance between two bins on a genome of genomeBins bins.
function circDist(a, b, genomeBins) {
    d = abs(a - b) % genomeBins;
    return minOf(d, genomeBins - d);
}

// Closes the window with the given title, if it's currently open. Loops
// rather than a single if-check since ImageJ does not enforce unique window
// titles, so more than one window can share the same title.
function closeIfOpen(title) {
    while (isOpen(title)) {
        selectWindow(title);
        close();
    }
}

// Adds two copies of `title`, shifted by +size and -size along `axis`
// ("x" or "y"), back into `title` -- makes the signal wrap periodically into
// the padded margin created in step 2, instead of fading to black at the edges.
function wrapAndAdd(title, size, negSize, axis) {
    shiftedPlus  = title + "_" + axis + "plus";
    shiftedMinus = title + "_" + axis + "minus";

    selectWindow(title);
    run("Duplicate...", "title=" + shiftedPlus);
    if (axis == "y")
        run("Translate...", "x=0 y=" + size + " interpolation=None");
    else
        run("Translate...", "x=" + size + " y=0 interpolation=None");

    selectWindow(title);
    run("Duplicate...", "title=" + shiftedMinus);
    if (axis == "y")
        run("Translate...", "x=0 y=" + negSize + " interpolation=None");
    else
        run("Translate...", "x=" + negSize + " y=0 interpolation=None");

    imageCalculator("Add", shiftedPlus, shiftedMinus);
    close(shiftedMinus);
    imageCalculator("Add", title, shiftedPlus);
    close(shiftedPlus);
}

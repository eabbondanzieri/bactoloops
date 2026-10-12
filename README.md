# Bactoloops

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23296823.svg)](https://doi.org/10.5281/zenodo.23296823)

Loop calling for Hi-C maps of circular bacterial chromosomes.

Published loop callers (HiCCUPS, cooltools, Chromosight, Mustache, SIP, Peakachu, ...) were built for linear eukaryotic
chromosomes. Bactoloops treats the bacterial chromosome as circular at every step:

- contact maps are padded periodically, so loops near the ends of the map are handled like any other;
- expected contacts, loop confidence and loop strength are computed from the circular genomic distance;
- the calls do not depend on where the circular map starts (rotating the map gives the same loops, confidence and scores);
- dense short-range contacts along the diagonal (chromosomal domains) are handled separately, with a stricter
  HiCCUPS-style donut test for loops close to the diagonal;
- each loop gets a confidence (its percentile among all pixels at the same circular distance) and an integrated excess
  (contacts above the local background, roughly proportional to how often the loop forms).

Bactoloops comes in two equivalent versions:

| | |
|---|---|
| `macro/Bactoloops.ijm` | ImageJ/Fiji macro (the reference implementation used in the paper) |
| `python/bactoloops.py` | Python implementation, no ImageJ needed; reproduces the macro's calls |

## Installation

**Fiji macro:** install [Fiji](https://fiji.sc); no plugins are needed. Open `macro/Bactoloops.ijm` in Fiji (Plugins >
Macros > Run...).

**Python:**

```bash
conda env create -f python/environment.yml
conda activate bactoloops
```

or `pip install -r python/requirements.txt`. The first run compiles part of the code with numba (a few seconds); a
461 x 461 map then takes about 20 seconds.

## Input

A square, ICE-normalized contact map saved as a 32-bit TIFF, one pixel per bin. The `tools/` folder converts the usual
Hi-C outputs:

```bash
# HiC-Pro valid pairs -> 10-kb ICE matrix (cis pairs, duplicates removed, pairs < 5 kb dropped)
python tools/validpairs_to_iced_matrix.py sample.validPairs sample_10000_iced.matrix --bin-size 10000 --genome-length 4607088

# sparse ICE matrix (HiC-Pro *_iced.matrix) -> TIFF, optionally rotated so that a chosen bin becomes the first
python tools/hic_circular_recenter_to_tiff.py --matrix sample_10000_iced.matrix --bin-size 10000 --out sample.tif --recenter-bin 368
```

The recentring only changes how the map is displayed. Give Bactoloops the same `--recenter-bin` (0 if not recentred)
so that it reports loops in genome coordinates.

## Running

**Fiji:** open the TIFF, run the macro, and fill in the dialog (output folder, identifier, recenter bin, bin size).
Batch use without the dialog is described in `docs/macro_details.md`.

**Python:**

```bash
python python/bactoloops.py sample.tif --recenter-bin 368 --bin-size 10000 --out results/ --id sample
```

## Output

Each file ends in the identifier you give:

| file | contents |
|---|---|
| `loop_calls_<id>.csv` | one row per loop: position (map bins and genome coordinates), confidence (high: > 99th percentile of pixels at the same distance; medium: > 95th), log2 observed/expected, integrated excess, and donut-test ratios |
| `integrated_loop_density_<id>.csv` | loop density summed along the genome (peaks mark loop anchors) |
| `Loop_density*.tif`, `Loop_density_smoothed_loops_<id>.png` | the processed map, with the called loops circled |
| `HiC_APA_<id>.tif` / `.png` | aggregate peak analysis of the called loops |
| `log_<id>.txt` | the run log, including every parameter |

The tunable parameters are listed at the top of the macro and of `bactoloops.py` (same names and defaults). They were
calibrated on 10-kb maps; `docs/macro_details.md` explains each step and how to run at other resolutions.

## Example and test

`example/ecoli_wt_96h_10kb.tif` is the *E. coli* wild-type deep-stationary-phase (96 h) map from the paper (average of two
biological replicates, 10-kb bins, recentred at bin 368). `example/expected/` holds the macro's output for it (128
loops). To check that the Python version reproduces the macro on your machine:

```bash
python tests/test_example.py
```

## Citation

If you use Bactoloops, please cite the software (all versions: doi:[10.5281/zenodo.23296823](https://doi.org/10.5281/zenodo.23296823); v1.0.0: doi:[10.5281/zenodo.23296824](https://doi.org/10.5281/zenodo.23296824)) and the paper:

Way, L. E. *et al.* H-NS organizes genome-wide DNA looping and enhances gene silencing in starved *Escherichia coli*.
bioRxiv 2025.10.14.682352 (2025); doi:10.1101/2025.10.14.682352.

The analysis scripts for the paper are in [eabbondanzieri/HNS-Deep-Stationary-Bridging-2026](https://github.com/eabbondanzieri/HNS-Deep-Stationary-Bridging-2026).

## License

MIT (see `LICENSE`).

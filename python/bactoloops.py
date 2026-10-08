#!/usr/bin/env python3
"""bactoloops.py -- Python version of HiC_PeakFindingAlgorithm_with_symmetrization.ijm.

Calls DNA loops on an ICE-normalized bacterial Hi-C contact map (a square float TIFF, one pixel
per bin), step for step as the ImageJ macro does, and writes the same outputs. The ImageJ
operations the macro relies on are reimplemented in ij_ops.py from the ImageJ 1.54 source; see
README_python.md for how closely this reproduces the macro and where it can differ.

Usage:
    python bactoloops.py map.tif --recenter-bin 368 --bin-size 10000 --out results/ --id wt_96h
Outputs (in --out, each ending in _<id>):
    loop_calls_<id>.csv                 one row per loop (same columns as the macro's table)
    Loop_density_<id>.tif               ridge-removed, noise-floor-subtracted map (16-bit)
    Loop_density_smoothed_<id>.tif      the smoothed map loops are called on (16-bit)
    HiC_APA_<id>.png / HiC_APA_<id>.tif aggregate peak analysis, log2(obs/exp) (PNG to view)
    integrated_loop_density_<id>.csv    genome_position_bp,ILD_au
    Loop_density_smoothed_loops_<id>.png the smoothed map with every loop circled
    log_<id>.txt                        the same summary lines the macro prints to its Log
Requires numpy, scipy, pandas, tifffile, numba, matplotlib (see environment.yml).
"""
import argparse
import os

import numpy as np
import pandas as pd
import tifffile

import ij_ops as ij

# ---- tunables: same names and defaults as the macro ----------------------------------------
targetRowSum = 5000000
padding = 200
medianRadius = 2
rollingBallRadius = 10
ridgeSmoothSigma = 3
minParticleSize = 2000
maskBlurSigma = 1
finalBlurSigma = 1
ridgeThresholdMultiplier = 0.5
noiseFloorMultiplier = 0.5
displayMax = 3000
minSmallParticleArea = 4
maxSmallParticleArea = 100
upsampleFactor = 5
robustThresholdMultiplier = 1.5
refineRadius = 1
localBackgroundRadius = 50
localBackgroundPhaseAverage = True   # rotation-invariant local background (see ij_ops.periodic_background)
localBackgroundPhases = 4
integratedExcessHalf = 2
confidenceHighPercentile = 0.99
confidenceMediumPercentile = 0.95
keepLowConfidence = False   # False: "low" loops are removed right after they are scored (as in the macro)
ridgeSearch = True         # look for loops inside the ridge mask (local maxima + HiCCUPS-style donut test)
ridgeMinPercentile = 0.95  # ridge search: candidates must be at least this distance-matched percentile
minLoopDist = 3            # bins -- no loop closer than this to the diagonal
donutPeakWidth = 2         # HiCCUPS -p at 10 kb: half-width of the peak square left out of every background
donutWindow = 5            # HiCCUPS -w at 10 kb: half-width of the donut window
donutRatioAll = 1.5        # obs / local expected must reach this for donut, lower-left, horizontal and vertical
donutRatioOne = 1.75       # ... and this for the donut or the lower-left (HiCCUPS thresholds, Rao et al. 2014)
apaSize = 9
loopCircleRadius = 6     # px (bins), as in the macro

LOG = []


def log(msg):
    LOG.append(msg)
    print(msg)


def circ_dist(a, b, n):
    d = abs(a - b) % n
    return min(d, n - d)


def undo_recenter(b, rc, n):
    return b if rc == 0 else (rc - 1 + b) % n


def undo_recenter_signed(b, rc, n):
    if rc == 0:
        return b
    raw = (rc - 1) + b
    g = raw % n
    if raw < n and rc > 1:
        g -= n
    return g


def ij_round(x):
    """ImageJ's round(): floor(x + 0.5) (Python's round() rounds halves to even)."""
    return int(np.floor(x + 0.5))


# ---- 7b ---------------------------------------------------------------------------------------
def resolve_symmetric_loops(P, n, rc, bin_bp):
    """resolveSymmetricLoops(): collapse the mirror-image detections of each loop (one on each
    side of the diagonal) into one row. P: particle rows with X, Y (centroid, image bins),
    Area, Mean, id."""
    if len(P) == 0:
        return pd.DataFrame()
    bx = np.floor(P.X.values).astype(int)
    by = np.floor(P.Y.values).astype(int)
    lo, hi = np.minimum(bx, by), np.maximum(bx, by)
    side = np.where(by < bx, 1, np.where(by > bx, -1, 0))
    clo, chi = np.minimum(P.X.values, P.Y.values), np.maximum(P.X.values, P.Y.values)
    key = lo * 100000 + hi
    order = np.argsort(key, kind="stable")
    out, pend = [], []
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and key[order[j + 1]] == key[order[i]]:
            j += 1
        grp = order[i:j + 1]
        sides = set(side[grp])
        if len(grp) >= 2 and (0 in sides or (1 in sides and -1 in sides)):
            out.append(dict(bin1=lo[grp[0]], bin2=hi[grp[0]], area=P.Area.values[grp].mean(),
                            mean=P.Mean.values[grp].mean(), ndet=len(grp), single=False,
                            ids=",".join(str(P.id.values[g]) for g in grp),
                            c1=clo[grp].mean(), c2=chi[grp].mean()))
        else:
            pend.extend(grp)
        i = j + 1
    used = [False] * len(pend)
    for a in range(len(pend)):
        if used[a]:
            continue
        for b in range(a + 1, len(pend)):
            if used[b]:
                continue
            pa, pb = pend[a], pend[b]
            if side[pa] != 0 and side[pb] != 0 and side[pa] == side[pb]:
                continue
            same_lo = lo[pa] == lo[pb] and abs(hi[pa] - hi[pb]) == 1
            same_hi = hi[pa] == hi[pb] and abs(lo[pa] - lo[pb]) == 1
            if same_lo or same_hi:
                out.append(dict(bin1=(lo[pa] + lo[pb]) / 2, bin2=(hi[pa] + hi[pb]) / 2,
                                area=(P.Area.values[pa] + P.Area.values[pb]) / 2,
                                mean=(P.Mean.values[pa] + P.Mean.values[pb]) / 2, ndet=2, single=False,
                                ids=f"{P.id.values[pa]},{P.id.values[pb]}",
                                c1=(clo[pa] + clo[pb]) / 2, c2=(chi[pa] + chi[pb]) / 2))
                used[a] = used[b] = True
                break
    for a, pa in enumerate(pend):
        if not used[a]:
            out.append(dict(bin1=lo[pa], bin2=hi[pa], area=P.Area.values[pa], mean=P.Mean.values[pa],
                            ndet=1, single=True, ids=str(P.id.values[pa]), c1=clo[pa], c2=chi[pa]))
    log(f"resolveSymmetricLoops: {len(P)} raw detections -> {len(out)} unique loops "
        f"({sum(o['single'] for o in out)} single-sided).")
    rows = []
    for o in out:
        g1, g2 = undo_recenter(o["bin1"], rc, n), undo_recenter(o["bin2"], rc, n)
        rows.append({"bin1": o["bin1"], "bin2": o["bin2"], "bp1": o["bin1"] * bin_bp, "bp2": o["bin2"] * bin_bp,
                     "genome_bin1": g1, "genome_bin2": g2, "genome_bp1": g1 * bin_bp, "genome_bp2": g2 * bin_bp,
                     "n_detections": o["ndet"],
                     "note": "TRUE single-sided detection" if o["single"] else "mirror (clean or rounding-boundary)",
                     "outline_ids": o["ids"], "mean_area": o["area"], "mean_intensity": o["mean"],
                     "loop_strength": o["area"] * o["mean"], "centroid1": o["c1"], "centroid2": o["c2"]})
    return pd.DataFrame(rows)


# ---- 7c helpers ----------------------------------------------------------------------------
def refine_pick(b1, b2, arr, n, radius, min_dist):
    """Offset (dy, dx) within `radius` of (b1, b2) of the highest value of arr; ties keep the
    unmoved pixel; moves closer than min_dist to the diagonal are skipped."""
    best = arr[b1 % n, b2 % n]
    s = (0, 0)
    for dy in range(-radius, radius + 1):
        for dx in range(-radius, radius + 1):
            if circ_dist(b1 + dy, b2 + dx, n) < min_dist:
                continue
            v = arr[(b1 + dy) % n, (b2 + dx) % n]
            if v > best:
                best, s = v, (dy, dx)
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tif")
    ap.add_argument("--recenter-bin", type=int, default=368,
                    help="the --recenter-bin used to build the TIFF (0 = not recentered)")
    ap.add_argument("--bin-size", type=int, default=10000, help="bp per pixel")
    ap.add_argument("--out", default=".")
    ap.add_argument("--id", default="")
    ap.add_argument("--no-ridge-search", action="store_true", help="skip the search inside the ridge mask")
    ap.add_argument("--dump-intermediates", action="store_true",
                    help="also save every intermediate image and candidate table (<out>/intermediates_<id>.npz, *_particles/_candidates/_scored .csv), e.g. for an explanatory figure")
    args = ap.parse_args()
    global ridgeSearch
    if args.no_ridge_search:
        ridgeSearch = False
    rc, bin_bp, ident = args.recenter_bin, args.bin_size, args.id
    os.makedirs(args.out, exist_ok=True)
    out = lambda name: os.path.join(args.out, f"{name}_{ident}")

    # ---- 0/1: load, scale, smooth the diagonal ----
    img = tifffile.imread(args.tif).astype(np.float32)
    n = img.shape[0]
    if not (0 <= rc <= n):
        raise SystemExit(f"Recenter bin must be 0 or between 1 and {n}; got {rc}.")
    row_sums = img.sum(axis=1, dtype=np.float64)
    row_sum = row_sums[np.argsort(row_sums, kind="stable")[n // 2]]
    scale = targetRowSum / row_sum
    log(f"HiC Loop Finding: row-sum = {row_sum}, scaleFactor = {scale}")
    img = (img * np.float32(scale)).astype(np.float32)
    log(f"HiC Loop Finding: {os.path.basename(args.tif)}, {n} bins, recenter bin = {rc}, bin size = {bin_bp} bp")
    if bin_bp != 10000:
        log("NOTE: the tunables that are set in pixels/bins were calibrated on 10 kb maps -- see the README.")
    # each diagonal pixel -> median of its 3x3 neighbourhood, wrapping around the circular map
    k = np.arange(n)
    nb = np.stack([img[(k + dy) % n, (k + dx) % n] for dy in (-1, 0, 1) for dx in (-1, 0, 1)])
    img[k, k] = np.sort(nb, axis=0)[4]

    # ---- 2: pad periodically ----
    ph = padding // 2
    main_ = np.pad(img, ph, mode="wrap").astype(np.float32)

    # ---- 3: background subtraction ----
    bg = ij.rolling_ball_background(ij.median(main_, medianRadius), rollingBallRadius)
    rb = (main_ - bg).astype(np.float32)

    # ---- 4: shot noise and the ridge mask ----
    # measured on the unpadded map only (the padding repeats parts of the map), as in the macro
    q1, med, q3 = ij.histogram_quantiles(rb[ph:ph + n, ph:ph + n], [0.25, 0.5, 0.75])
    sd = (q3 - q1) / 1.349
    log(f"Shot-noise estimate (32-bit HiC_map_RB, pre-ridge-removal) -- median={med}, robustSD={sd}")
    sm = ij.gaussian(ij.median(rb, medianRadius), ridgeSmoothSigma)
    thr = med + ridgeThresholdMultiplier * sd
    lab, keep = ij.particles(sm >= thr, min_area=minParticleSize)
    ridge = np.isin(lab, keep).astype(np.float32) * 255
    ridge = ij.gaussian(ridge, maskBlurSigma, bit_depth=8)
    mask = (ridge < 1).astype(np.float32)             # 0 inside the ridge, 1 outside (Max 1, x256, Invert, /256)
    loop_density = (rb * mask).astype(np.float32)

    # ---- 5: noise floor, 16-bit ----
    offset = med + noiseFloorMultiplier * sd
    log(f"Step 5 shot-noise floor: median={med}, robustSD={sd}, noiseFloorMultiplier={noiseFloorMultiplier} -> subtracting {offset}")
    loop_density = np.clip(np.floor(loop_density - offset + 0.5), 0, 65535).astype(np.float32)

    # ---- 6: Loop_density_smoothed: loops are called on the masked version; the saved/displayed
    # version shows the data everywhere (same steps from rb, before the ridge mask)
    smoothed_pad = ij.gaussian(loop_density, finalBlurSigma, bit_depth=16)   # padded: particles at the map edge are found whole
    smoothed = smoothed_pad[ph:ph + n, ph:ph + n]
    display = np.clip(np.floor(rb - offset + 0.5), 0, 65535).astype(np.float32)
    display = ij.gaussian(display, finalBlurSigma, bit_depth=16)[ph:ph + n, ph:ph + n]

    # ---- 7: particles on the upsampled, padded map; keep those centred inside the map ----
    up = ij.resize_bilinear_16bit(smoothed_pad, upsampleFactor)
    thr7 = 0 + (robustThresholdMultiplier - noiseFloorMultiplier) * sd
    binary = ij.watershed_split(up >= thr7)
    lab, keep = ij.particles(binary, minSmallParticleArea, maxSmallParticleArea, area_scale=1 / upsampleFactor ** 2)
    rows = []
    for k, L in enumerate(keep):
        yy, xx = np.nonzero(lab == L)
        X, Y = (xx.mean() + 0.5) / upsampleFactor, (yy.mean() + 0.5) / upsampleFactor
        if ph <= X < ph + n and ph <= Y < ph + n:        # the central copy; copies in the padding are dropped
            rows.append(dict(id=len(rows), X=X - ph, Y=Y - ph, Area=len(xx) / upsampleFactor ** 2, Mean=up[yy, xx].mean()))
    log(f"Step 7: {len(keep)} particles on the padded map, {len(rows)} with their centroid inside the map.")
    P = pd.DataFrame(rows)

    # ---- 7b ----
    T = resolve_symmetric_loops(P, n, rc, bin_bp)
    T_candidates = T.copy()                      # before refinement (for --dump-intermediates)

    # ---- 7c: obs/exp, local background, refinement, confidence, ridge search, excess ----
    nd = n // 2 + 1
    core = main_[ph:ph + n, ph:ph + n].astype(np.float64)
    i_, j_ = np.indices((n, n))
    dist = np.minimum(abs(i_ - j_), n - abs(i_ - j_))
    upper = (j_ >= i_) & (core > 0)
    esum = np.bincount(dist[upper], core[upper], minlength=nd)
    ecnt = np.bincount(dist[upper], minlength=nd)
    expected = np.where(ecnt > 0, esum / np.maximum(ecnt, 1), 0.0)
    has = core > 0
    nonlocal_rows = np.where(dist >= 2, np.where(has, core, 0), 0).sum(1)
    nonlocal_med = np.median(nonlocal_rows[has.any(1)])
    # log2(obs/exp) on the padded canvas (periodic), missing -> 0, 3x3 average
    pi, pj = np.indices(main_.shape)
    bi, bj = (pi - ph) % n, (pj - ph) % n
    pdist = np.minimum(abs(bi - bj), n - abs(bi - bj))
    E = np.where(expected[pdist] > 0, expected[pdist], 1.0).astype(np.float32)
    oe = (main_ / E).astype(np.float32)
    oe = np.where(~(oe > 0) | (oe > 1e30), 1.0, oe).astype(np.float32)
    log1 = (np.log(oe) * np.float32(1 / np.log(2))).astype(np.float32)
    kern = np.ones((3, 3), np.float32) / 9
    from scipy import ndimage as ndi
    avg_raw = ndi.correlate(log1, kern, mode="nearest").astype(np.float32)
    if localBackgroundRadius <= 0:
        avg = avg_raw
    elif localBackgroundPhaseAverage:
        avg = (avg_raw - ij.periodic_background(avg_raw, n, ph, localBackgroundRadius, localBackgroundPhases)).astype(np.float32)
    else:
        avg = ij.subtract_background(avg_raw, localBackgroundRadius)
    C = lambda a: a[ph:ph + n, ph:ph + n].astype(np.float64)
    rawA, avgRawA, avgA = C(log1), C(avg_raw), C(avg)
    bg_removed = avgRawA - avgA
    refA = rawA - bg_removed
    excA = np.where(core > 0, core - expected[dist] * 2.0 ** bg_removed, 0.0)

    if len(T):
        T["refine_shift1"] = 0
        T["refine_shift2"] = 0
    n_moved = n_merged = 0
    if refineRadius > 0 and len(T):
        for r in T.index:
            b1, b2 = ij_round(T.at[r, "bin1"]), ij_round(T.at[r, "bin2"])
            s1, s2 = refine_pick(b1, b2, refA, n, refineRadius, minLoopDist)
            n1, n2 = (b1 + s1) % n, (b2 + s2) % n
            n_moved += (s1, s2) != (0, 0)
            T.at[r, "bin1"], T.at[r, "bin2"] = min(n1, n2), max(n1, n2)
            T.at[r, "refine_shift1"], T.at[r, "refine_shift2"] = s1, s2
        keep_rows = []
        for r in T.index:                       # merge loops that land on the same pixel (keep the first)
            dup = [k for k in keep_rows if T.at[k, "bin1"] == T.at[r, "bin1"] and T.at[k, "bin2"] == T.at[r, "bin2"]]
            if dup:
                k = dup[0]
                T.at[k, "outline_ids"] = f"{T.at[k, 'outline_ids']},{T.at[r, 'outline_ids']}"
                T.at[k, "n_detections"] += T.at[r, "n_detections"]
                n_merged += 1
            else:
                keep_rows.append(r)
        T = T.loc[keep_rows].reset_index(drop=True)
        for r in T.index:
            b1, b2 = T.at[r, "bin1"], T.at[r, "bin2"]
            g1, g2 = undo_recenter(int(b1), rc, n), undo_recenter(int(b2), rc, n)
            T.loc[r, ["bp1", "bp2", "genome_bin1", "genome_bin2", "genome_bp1", "genome_bp2"]] = \
                [b1 * bin_bp, b2 * bin_bp, g1, g2, g1 * bin_bp, g2 * bin_bp]
        log(f"Step 7c refinement: {n_moved} of {len(T) + n_merged} loops moved to the single-pixel log2(obs/exp) "
            f"maximum within {refineRadius} bin(s); {n_merged} merged after landing on the same pixel.")

    def null_values(d):
        """values the loop is ranked against: every pixel at distance d, or (with refinement) the
        3x3-averaged value at each such pixel's own refinement pick"""
        ks = np.arange(n)
        if refineRadius <= 0:
            return avgA[ks, (ks + d) % n]
        vals = np.empty(n)
        for k in ks:
            kc = (k + d) % n
            s1, s2 = refine_pick(k, kc, refA, n, refineRadius, minLoopDist)
            vals[k] = avgA[(k + s1) % n, (kc + s2) % n]
        return vals

    def confidence(p):
        return "high" if p >= confidenceHighPercentile else "medium" if p >= confidenceMediumPercentile else "low"

    counts = {"high": 0, "medium": 0, "low": 0}
    for r in T.index:
        b1, b2 = ij_round(T.at[r, "bin1"]), ij_round(T.at[r, "bin2"])
        val = avgA[b1, b2]
        pct = np.mean(null_values(circ_dist(b1, b2, n)) < val)
        T.at[r, "log2_obs_exp"], T.at[r, "oe_percentile"] = val, pct
        T.at[r, "confidence"] = confidence(pct)
        counts[T.at[r, "confidence"]] += 1
    log(f"Step 7c confidence: {counts['high']} high, {counts['medium']} medium, {counts['low']} low.")
    T_scored = T.copy()                          # every candidate with its confidence, low included
    if not keepLowConfidence and counts["low"] > 0:
        T = T[T.confidence != "low"].reset_index(drop=True)
        log(f"Step 7c: removed the {counts['low']} low-confidence loop(s) (keepLowConfidence = False); {len(T)} remain.")

    ridge_core = mask[ph:ph + n, ph:ph + n]
    null_cache = {}

    def pct_at(b1, b2):
        d = circ_dist(b1, b2, n)
        if d not in null_cache:
            null_cache[d] = null_values(d)
        return float(np.mean(null_cache[d] < avgA[b1, b2]))

    # HiCCUPS-style local backgrounds (observed / locally expected) for every pixel we test.
    obs = np.where(core > 0, core, np.nan)
    expd = expected[dist]
    def donut_ratios(i, j):
        """oriented along the short arc (j at most n/2 bins downstream of i), so the lower-left
        kernel faces the diagonal whatever the map's start; a pixel exactly half the genome away
        gets the mean of both orientations"""
        fwd = (j - i) % n
        if 2 * fwd == n:
            a, b = donut_ratios_oriented(i, j), donut_ratios_oriented(j, i)
            return {k: (a[k] + b[k]) / 2 for k in a}
        return donut_ratios_oriented(j, i) if fwd > n - fwd else donut_ratios_oriented(i, j)

    def donut_ratios_oriented(i, j):
        p, w = donutPeakWidth, donutWindow
        rows = np.arange(i - w, i + w + 1); cols = np.arange(j - w, j + w + 1)
        di, dj = np.meshgrid(rows - i, cols - j, indexing="ij")
        O = obs[np.ix_(rows % n, cols % n)]; Ed = expd[np.ix_(rows % n, cols % n)]
        dd = np.minimum(abs((rows[:, None] % n) - (cols[None, :] % n)), n - abs((rows[:, None] % n) - (cols[None, :] % n)))
        ok = np.isfinite(O) & (dd >= 2)                      # leave out missing pixels and the diagonal
        peak = (abs(di) <= p) & (abs(dj) <= p)
        kern = {"donut": ~peak & (di != 0) & (dj != 0),
                "lowerleft": (di > 0) & (dj < 0) & ~peak,     # towards the diagonal (row > i, column < j)
                "horizontal": (abs(di) <= 1) & (abs(dj) > p),
                "vertical": (abs(dj) <= 1) & (abs(di) > p)}
        centre = obs[i % n, j % n]
        out = {}
        for k, m in kern.items():
            m = m & ok
            if m.sum() == 0 or not np.isfinite(centre):
                out[k] = np.nan
                continue
            local_exp = O[m].sum() / Ed[m].sum() * expd[i % n, j % n]
            out[k] = centre / local_exp
        return out

    def passes_donut(r):
        vals = [r["donut"], r["lowerleft"], r["horizontal"], r["vertical"]]
        return (all(np.isfinite(v) and v >= donutRatioAll for v in vals)
                and max(r["donut"], r["lowerleft"]) >= donutRatioOne)

    # ridge search 1 (default): every local maximum inside the ridge mask, tested against its
    # distance-matched null and then the donut backgrounds
    if ridgeSearch:
        m1 = np.array([ij_round(v) for v in T.bin1]) if len(T) else np.array([], int)
        m2 = np.array([ij_round(v) for v in T.bin2]) if len(T) else np.array([], int)
        cands, n_max, n_pct = [], 0, 0
        for i in range(n):
            for d in range(minLoopDist, n // 2 + 1):
                j = (i + d) % n
                if ridge_core[i, j] >= 0.5:
                    continue                                  # outside the ridge: the main pass covers it
                v = refA[i, j]
                nb = refA[np.ix_(np.arange(i - 1, i + 2) % n, np.arange(j - 1, j + 2) % n)]
                if v < nb.max():
                    continue                                  # not a local maximum
                n_max += 1
                b1, b2 = min(i, j), max(i, j)
                if any(abs(b1 - x) <= 3 and abs(b2 - y) <= 3 for x, y in zip(m1, m2)):
                    continue                                  # already a main-pass loop
                pct = pct_at(i, j)
                if pct < ridgeMinPercentile:
                    continue
                n_pct += 1
                r = donut_ratios(i, j)
                if passes_donut(r):
                    cands.append((avgA[i, j], b1, b2, pct, r))
        added = []
        for val, b1, b2, pct, r in sorted(cands, key=lambda t: t[0], reverse=True):
            if any(abs(b1 - x) <= 3 and abs(b2 - y) <= 3 for x, y, *_ in added):
                continue
            added.append((b1, b2, val, pct, r))
        new = []
        for b1, b2, val, pct, r in added:
            g1, g2 = undo_recenter(b1, rc, n), undo_recenter(b2, rc, n)
            new.append({"bin1": b1, "bin2": b2, "bp1": b1 * bin_bp, "bp2": b2 * bin_bp, "genome_bin1": g1,
                        "genome_bin2": g2, "genome_bp1": g1 * bin_bp, "genome_bp2": g2 * bin_bp, "n_detections": 0,
                        "note": "ridge region (donut test)", "outline_ids": "", "mean_area": np.nan,
                        "mean_intensity": np.nan, "loop_strength": np.nan, "centroid1": np.nan, "centroid2": np.nan,
                        "refine_shift1": 0, "refine_shift2": 0, "log2_obs_exp": val, "oe_percentile": pct,
                        "confidence": confidence(pct)})
        if new:
            T = pd.concat([T, pd.DataFrame(new)], ignore_index=True)
        log(f"Step 7c ridge search: {n_max} local maxima inside the ridge mask (>= {minLoopDist} bins from the diagonal); "
            f"{n_pct} at or above the {100 * ridgeMinPercentile:.0f}th distance-matched percentile; "
            f"{len(cands)} passed the donut test; {len(added)} added after merging calls within 3 bins.")

    if len(T):                       # donut ratios for every loop (a check; only ridge calls must pass)
        R = [donut_ratios(ij_round(a), ij_round(b)) for a, b in zip(T.bin1, T.bin2)]
        for k in ["donut", "lowerleft", "horizontal", "vertical"]:
            T[f"{k}_ratio"] = [r[k] for r in R]
    if len(T):
        h = integratedExcessHalf
        T["integrated_excess"] = [100 * excA[np.ix_(np.arange(ij_round(a) - h, ij_round(a) + h + 1) % n,
                                                     np.arange(ij_round(b) - h, ij_round(b) + h + 1) % n)].sum() / nonlocal_med
                                  for a, b in zip(T.bin1, T.bin2)]
    log(f"Step 7c: integrated_excess summed over {2 * integratedExcessHalf + 1} x {2 * integratedExcessHalf + 1} bins, "
        f"as % of the median bin's non-local contacts ({nonlocal_med}).")

    # ---- 7e: APA ----
    ah = (apaSize - 1) // 2
    apa_sum = np.zeros((apaSize, apaSize)); apa_n = np.zeros((apaSize, apaSize))
    for a, b in zip(T.bin1 if len(T) else [], T.bin2 if len(T) else []):
        r0, c0 = ij_round(a), ij_round(b)
        for dy in range(-ah, ah + 1):
            for dx in range(-ah, ah + 1):
                v = main_[r0 + dy + ph, c0 + dx + ph]
                d = circ_dist(c0 + dx, r0 + dy, n)
                if v > 0 and expected[d] > 0:
                    apa_sum[dy + ah, dx + ah] += v / expected[d]; apa_n[dy + ah, dx + ah] += 1
    apa = np.where(apa_n > 0, np.log2(np.where(apa_n > 0, apa_sum / np.maximum(apa_n, 1), 1)), 0).astype(np.float32)
    log(f"Step 7e APA: averaged {len(T)} loop(s); centre pixel log2(obs/exp) = {apa[ah, ah]}")

    # ---- 8: integrated loop density ----
    ld_core = loop_density[ph:ph + n, ph:ph + n]
    profile = ld_core.mean(axis=0)
    xbp = [undo_recenter_signed(i, rc, n) * bin_bp for i in range(n)]

    # ---- 9: save ----
    tifffile.imwrite(out("Loop_density") + ".tif", ld_core.astype(np.uint16))
    tifffile.imwrite(out("Loop_density_smoothed") + ".tif", display.astype(np.uint16))
    tifffile.imwrite(out("HiC_APA") + ".tif", apa)
    cols = ["bin1", "bin2", "bp1", "bp2", "genome_bin1", "genome_bin2", "genome_bp1", "genome_bp2", "n_detections",
            "note", "outline_ids", "mean_area", "mean_intensity", "loop_strength", "centroid1", "centroid2",
            "refine_shift1", "refine_shift2", "log2_obs_exp", "oe_percentile", "confidence", "integrated_excess",
            "donut_ratio", "lowerleft_ratio", "horizontal_ratio", "vertical_ratio"]
    (T[cols] if len(T) else pd.DataFrame(columns=cols)).to_csv(out("loop_calls") + ".csv", index=False)
    if args.dump_intermediates:
        np.savez_compressed(out("intermediates") + ".npz", n=n, pad_half=ph, upsample=upsampleFactor,
                            scaled_padded=main_, rolling_ball_bg=bg, bg_subtracted=rb, ridge_smoothed=sm,
                            ridge_keep_mask=mask, loop_density=loop_density, loop_density_smoothed_padded=smoothed_pad,
                            display=display, particles_binary=binary, shot_noise_median=med, shot_noise_sd=sd,
                            particle_threshold=thr7, expected=expected, log2oe_raw=rawA, log2oe_3x3=avgRawA,
                            log2oe_corrected=avgA, refine_map=refA, excess=excA, nonlocal_median=nonlocal_med)
        P.to_csv(out("particles") + ".csv", index=False)
        T_candidates.to_csv(out("candidates") + ".csv", index=False)
        T_scored.to_csv(out("scored") + ".csv", index=False)
    pd.DataFrame({"genome_position_bp": xbp, "ILD_au": profile}).to_csv(out("integrated_loop_density") + ".csv", index=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    lim = max(abs(apa.min()), abs(apa.max())) or 1
    rdbu = LinearSegmentedColormap.from_list("rdbu", ["#2166ac", "#ffffff", "#b2182b"])
    plt.imsave(out("HiC_APA") + ".png", np.kron(apa, np.ones((40, 40))), cmap=rdbu, vmin=-lim, vmax=lim)
    fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
    ax.imshow(display, cmap="gray_r", vmin=0, vmax=displayMax)
    ax.contour((ridge_core < 0.5).astype(float), levels=[0.5], colors=["#00c8e0"], linewidths=1.0)
    for _, r in T.iterrows():
        low = r.confidence == "low"
        ax.add_patch(plt.Circle((r.bin2, r.bin1), loopCircleRadius, fill=False, lw=0.8,
                                color="#99b8ff" if low else "blue", ls="--" if low else "-"))
        if str(r.note).startswith("ridge"):
            ax.plot(r.bin2, r.bin1, "+", color="green", ms=5, mew=0.8)
    ax.axis("off")
    fig.savefig(out("Loop_density_smoothed_loops") + ".png", bbox_inches="tight")
    with open(out("log") + ".txt", "w") as fh:
        fh.write("\n".join(LOG) + "\n")
    log(f"Saved outputs to {args.out} with identifier '{ident}'.")


if __name__ == "__main__":
    main()

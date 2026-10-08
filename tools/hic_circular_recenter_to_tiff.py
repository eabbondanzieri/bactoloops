#!/usr/bin/env python3
"""
hic_circular_recenter_to_tiff.py

Python port of hic_circular_recenter_to_tiff.R: bare-bones converter for a
single circular bacterial chromosome. Reads one sparse Hi-C contact matrix
file (bin1  bin2  value) at a known resolution and writes it out as a
single-channel 32-bit floating-point TIFF -- raw values, no rescaling,
clipping, normalization, or diagonal masking.

The only optional post-processing step is a circular "recenter": rotate the
matrix so a chosen bin becomes the new bin 1, so a feature of interest
(e.g. the replication terminus) lands in the middle of the plot instead of
at the corners. Leave --recenter-bin unset to skip this and keep the
matrix in its original bin order.

USAGE
-----
    # Original bin order, no recentering:
    python3 hic_circular_recenter_to_tiff.py --matrix genome_10000_iced.matrix \\
        --bin-size 10000 --out contacts.tif

    # Recentered so original bin 368 becomes the new bin 1:
    python3 hic_circular_recenter_to_tiff.py --matrix genome_10000_iced.matrix \\
        --bin-size 10000 --out contacts_recentered.tif --recenter-bin 368

    Other optional flags: [--log1p] [--id-offset 0] [--id-mode auto|bp|index]
                          [--max-pixels 4e7]

Requires: numpy, pandas, tifffile (pip install numpy pandas tifffile)
"""

import argparse
import sys

import numpy as np
import pandas as pd
import tifffile


def build_and_convert(matrix_file, bin_size, out_file,
                       id_mode="auto",
                       log1p_transform=False,
                       id_offset=0,
                       max_pixels=4e7,
                       recenter_bin=None):

    mat = pd.read_csv(matrix_file, sep=r"\s+", header=None,
                       names=["bin1", "bin2", "value"])

    # --- validate: catch malformed/non-numeric files with a clear message ---
    # (e.g. a stray non-matrix file matched by a glob pattern in a batch run)
    if mat.shape[1] != 3:
        sys.exit(f"{matrix_file}: expected 3 whitespace-separated columns "
                  f"(bin1 bin2 value), got {mat.shape[1]}")
    for col in ("bin1", "bin2", "value"):
        mat[col] = pd.to_numeric(mat[col], errors="coerce")
    bad_rows = mat[["bin1", "bin2", "value"]].isna().any(axis=1)
    if bad_rows.any():
        first_bad = mat.loc[bad_rows].iloc[0]
        sys.exit(f"{matrix_file}: {int(bad_rows.sum())} row(s) aren't numeric "
                  f"bin1/bin2/value triplets -- is this really a sparse Hi-C "
                  f"matrix file? First bad row (0-indexed row {bad_rows.idxmax()}): "
                  f"{first_bad.to_dict()}")

    # --- work out bin1/bin2 -> matrix row/col index -------------------------
    if id_mode == "auto":
        divisible = bin_size > 1 and \
            (mat["bin1"] % bin_size == 0).all() and \
            (mat["bin2"] % bin_size == 0).all()
        id_mode = "bp" if divisible else "index"
        print(f"Auto-detected --id-mode = '{id_mode}'")

    if id_mode == "bp":
        # bin1/bin2 are genomic bp start coordinates, multiples of bin_size
        i_idx = np.round(mat["bin1"].to_numpy() / bin_size).astype(np.int64) + 1 + id_offset
        j_idx = np.round(mat["bin2"].to_numpy() / bin_size).astype(np.int64) + 1 + id_offset
    elif id_mode == "index":
        # bin1/bin2 are already small sequential bin indices
        i_idx = mat["bin1"].to_numpy().astype(np.int64) + id_offset
        j_idx = mat["bin2"].to_numpy().astype(np.int64) + id_offset
    else:
        sys.exit(f"--id-mode must be one of 'auto', 'bp', or 'index', got: {id_mode}")

    if (i_idx < 1).any() or (j_idx < 1).any():
        sys.exit(
            f"Computed bin indices go below 1 (min i={i_idx.min()}, "
            f"min j={j_idx.min()}). Adjust --id-offset (try 0 or 1) or "
            f"check --id-mode."
        )

    n = int(max(i_idx.max(), j_idx.max()))
    if n ** 2 > max_pixels:
        sys.exit(
            f"Matrix would be {n} x {n} bins ({n**2/1e6:.1f}M pixels), over "
            f"--max-pixels ({max_pixels/1e6:.1f}M). Raise --max-pixels if you "
            f"really have enough RAM."
        )
    print(f"Building a {n} x {n} contact matrix "
          f"(genome length ~{n * bin_size / 1e6:.3f} Mb at {bin_size} bp bins)")

    # --- normalize to one triangle, collapsing (i,j)/(j,i) duplicates ------
    # Matrix files vary: some store only the upper triangle, some only the
    # lower, some a mix of both (a data-quality quirk), and some store BOTH
    # (i,j) and (j,i) for every pair. Folding everything onto lo<=hi and
    # averaging any duplicates handles all of these without assuming which
    # convention your file uses.
    lo = np.minimum(i_idx, j_idx)
    hi = np.maximum(i_idx, j_idx)
    values = mat["value"].to_numpy(dtype=np.float64)

    dt = pd.DataFrame({"lo": lo, "hi": hi, "value": values})
    grouped = dt.groupby(["lo", "hi"])["value"].agg(["mean", "count", "max", "min"])
    grouped = grouped.reset_index()
    grouped["spread"] = np.where(grouped["count"] > 1,
                                  grouped["max"] - grouped["min"], 0.0)

    n_collapsed = int((grouped["count"] > 1).sum())
    if n_collapsed > 0:
        max_rel_spread = (grouped["spread"] / np.maximum(np.abs(grouped["mean"]), 1e-12)).max()
        print(
            f"Note: {n_collapsed} bin pairs appeared more than once (both "
            f"directions and/or duplicate rows) and were averaged (max "
            f"relative difference among duplicates: {100 * max_rel_spread:.2g}%)."
        )
        if max_rel_spread > 0.05:
            print(
                "WARNING: Duplicate entries for the same bin pair disagree by "
                "more than 5% in places -- double check the matrix file isn't "
                "mixing two different runs/normalizations.", file=sys.stderr
            )

    # --- assemble the dense symmetric matrix -------------------------------
    m = np.zeros((n, n), dtype=np.float64)
    lo0 = grouped["lo"].to_numpy() - 1  # to 0-based
    hi0 = grouped["hi"].to_numpy() - 1
    vals = grouped["mean"].to_numpy()
    m[lo0, hi0] = vals
    m = m + m.T - np.diag(np.diag(m))  # mirror without doubling the diagonal

    if log1p_transform:
        m = np.log1p(m)

    # --- optional circular recenter (2x2 tile + crop) -----------------------
    if recenter_bin is not None:
        if recenter_bin < 1 or recenter_bin > n:
            sys.exit(f"--recenter-bin must be between 1 and {n}, got {recenter_bin}")
        g = np.block([[m, m], [m, m]])
        start = recenter_bin - 1  # to 0-based
        m = g[start:start + n, start:start + n]
        print(f"Recentered: circularly rotated the matrix so original bin "
              f"{recenter_bin} is now bin 1")

    print(f"Value range written to TIFF: [{np.nanmin(m):.6g}, {np.nanmax(m):.6g}]")

    # --- write as a true 32-bit float TIFF (no additional scaling) ---------
    tifffile.imwrite(out_file, m.astype(np.float32))
    print(f"Wrote 32-bit float TIFF: {out_file}")

    return m


def main():
    p = argparse.ArgumentParser(
        description="Convert a single-chromosome (circular bacterial genome) "
                    "Hi-C sparse matrix file to a 32-bit float TIFF, with an "
                    "optional circular recenter.")
    p.add_argument("--matrix", required=True, help="Sparse matrix file (bin1 bin2 value)")
    p.add_argument("--bin-size", type=int, required=True, help="Resolution in bp")
    p.add_argument("--out", required=True, help="Output TIFF path")
    p.add_argument("--id-mode", choices=["auto", "bp", "index"], default="auto")
    p.add_argument("--id-offset", type=int, default=0)
    p.add_argument("--max-pixels", type=float, default=4e7)
    p.add_argument("--log1p", action="store_true")
    p.add_argument("--recenter-bin", type=int, default=None)
    args = p.parse_args()

    build_and_convert(
        matrix_file=args.matrix,
        bin_size=args.bin_size,
        out_file=args.out,
        id_mode=args.id_mode,
        log1p_transform=args.log1p,
        id_offset=args.id_offset,
        max_pixels=args.max_pixels,
        recenter_bin=args.recenter_bin,
    )


if __name__ == "__main__":
    main()

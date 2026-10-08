"""Build a HiC-Pro-style 10-kb ICE matrix from a per-library validPairs file (as HiC-Pro would from allValidPairs).
Steps (HiC-Pro defaults used in the lab pipeline): cis pairs, duplicates removed (same chr/pos/strand for both ends), pairs closer than
MIN_CIS_DIST = 5 kb dropped, bins = floor(pos / binsize) (1-based in the output), symmetric raw matrix; ICE as in iced
(filter_low_counts_perc 0.02 on the number of non-zero contacts per bin, remove all-zero loci, max_iter 100, eps 0.1, l1: total counts
kept). Output: sparse upper-triangle 'bin1 bin2 value' (.matrix), the same format as HiC-Pro's *_iced.matrix.
Usage: python validpairs_to_iced_matrix.py in.validPairs out_10000_iced.matrix [--bin-size 10000] [--genome-length 4607088]
(--genome-length: length of the circular chromosome in bp; the default is the E. coli W3110 genome used in our paper, with the
1,398,953-1,438,196 deletion of the lab strain.)"""
import argparse, numpy as np, pandas as pd
ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
ap.add_argument("validpairs"); ap.add_argument("out")
ap.add_argument("bin_size_pos", nargs="?", type=int, help=argparse.SUPPRESS)
ap.add_argument("--bin-size", type=int, default=10000); ap.add_argument("--genome-length", type=int, default=4_607_088)
a = ap.parse_args(); vp, out = a.validpairs, a.out; BS = a.bin_size_pos or a.bin_size
G = a.genome_length; NB = -(-G // BS)
P = pd.read_csv(vp, sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 6], names=["c1", "p1", "s1", "c2", "p2", "s2"])
P = P[P.c1 == P.c2].drop_duplicates(["p1", "s1", "p2", "s2"])
P = P[(P.p2 - P.p1).abs() >= 5000]
i, j = (P.p1.values // BS), (P.p2.values // BS)
M = np.zeros((NB, NB)); np.add.at(M, (i, j), 1); np.add.at(M, (j, i), 1)
M[np.arange(NB), np.arange(NB)] /= 2                       # a pair inside one bin was added twice
print(f"pairs used: {len(P):,}")
X = M.copy(); total = np.triu(X).sum()
nz = (X != 0).sum(0)
bad = nz == 0
thr = np.percentile(nz[~bad], 2)                           # filter_low_counts_perc 0.02 (sparsity)
bad |= nz < thr
X[bad, :] = 0; X[:, bad] = 0
bias = np.ones(NB)
for it in range(100):
    s = X.sum(0); m = s[s != 0].mean(); s = s / m; s[s == 0] = 1
    X = X / s[:, None] / s[None, :]; bias *= s
    if np.abs(s - 1).max() < 0.1: break
X *= total / np.triu(X).sum()
print(f"ICE: {it+1} iterations, {int(bad.sum())} bins filtered")
r, c = np.triu_indices(NB); v = X[r, c]; k = v != 0
pd.DataFrame({"a": r[k] + 1, "b": c[k] + 1, "v": v[k]}).to_csv(out, sep="\t", header=False, index=False, float_format="%.6f")

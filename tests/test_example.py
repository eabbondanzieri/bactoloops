"""Check that the Python implementation reproduces the ImageJ macro on the example map.
Runs python/bactoloops.py on example/ecoli_wt_96h_10kb.tif and compares the loop calls with the macro's output in
example/expected/: the same loops at the same positions, the same confidence, and scores within the macro's 3-decimal output.
Usage: python tests/test_example.py        (or: pytest tests/)"""
import os, subprocess, sys, tempfile
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.path.dirname(HERE)


def test_example_matches_macro():
    with tempfile.TemporaryDirectory() as out:
        subprocess.run([sys.executable, os.path.join(REPO, "python", "bactoloops.py"), os.path.join(REPO, "example", "ecoli_wt_96h_10kb.tif"),
                        "--recenter-bin", "368", "--bin-size", "10000", "--out", out, "--id", "test"], check=True)
        py = pd.read_csv(os.path.join(out, "loop_calls_test.csv"))
    mac = pd.read_csv(os.path.join(REPO, "example", "expected", "loop_calls_macro_ecoli_wt_96h.csv"))
    key = ["genome_bin1", "genome_bin2"]
    assert len(py) == len(mac), (len(py), len(mac))
    j = mac.set_index(key).join(py.set_index(key), rsuffix="_py", how="inner")
    assert len(j) == len(mac), "loop positions differ"
    assert (j.confidence == j.confidence_py).all(), "confidence differs"
    for c in ["log2_obs_exp", "oe_percentile", "integrated_excess"]:
        assert (j[c] - j[c + "_py"]).abs().max() <= 1e-3, c
    print(f"OK: {len(py)} loops, identical positions and confidence; scores within 0.001 of the macro")


if __name__ == "__main__":
    test_example_matches_macro()

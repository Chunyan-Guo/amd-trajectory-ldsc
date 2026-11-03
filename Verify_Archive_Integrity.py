#!/usr/bin/env python3
"""Verify_Archive_Integrity.py - run from the repository root after extraction.

Checks:
  1. every file listed in SHA256SUMS.txt exists and matches its sha256;
  2. no unexpected extra files (beyond MANIFEST.csv, SHA256SUMS.txt, this script);
  3. required tree elements are present (scripts, config, all 20 result JSONs).
Prints PASSED only if every check succeeds.
"""
import hashlib, os, sys

ROOT = os.path.dirname(os.path.abspath(__file__))

def sha256(path, _buf=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(_buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def main():
    sums_file = os.path.join(ROOT, "SHA256SUMS.txt")
    if not os.path.exists(sums_file):
        sys.exit("FAILED: SHA256SUMS.txt missing")
    listed = []
    bad = []
    with open(sums_file) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            digest, rel = line.split("  ", 1)
            listed.append(rel)
            fp = os.path.join(ROOT, rel)
            if not os.path.exists(fp):
                bad.append(("missing", rel))
            elif sha256(fp) != digest:
                bad.append(("hash mismatch", rel))
    if bad:
        for kind, rel in bad:
            print("BAD:", kind, rel)
        sys.exit(f"FAILED: {len(bad)} file(s)")

    allowed = {"MANIFEST.csv", "SHA256SUMS.txt", "Verify_Archive_Integrity.py"}
    present = set()
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if d not in {"__pycache__", ".git"}]
        for fn in filenames:
            rel = os.path.relpath(os.path.join(dirpath, fn), ROOT).replace(os.sep, "/")
            present.add(rel)
    extra = sorted(present - set(listed) - allowed)
    if extra:
        for rel in extra:
            print("UNLISTED:", rel)
        sys.exit("FAILED: %d unlisted file(s) in tree" % len(extra))

    required = (
        ["README.md", "LICENSE", "requirements.txt", "config/endpoints.tsv"]
        + [f"scripts/{n}" for n in ["00_check_environment.py", "01_download.py",
                                    "02_munge.py", "03_rg.py", "common.py"]]
        + [f"analysis/{n}" for n in ["emit_numbers.py", "fig1.py", "fig2.py", "make_supp.py"]]
        + [f"results/qc_{t}.json" for t in
           ["mac", "amd", "dry", "wet", "htn", "str", "t2d", "ckd", "glom", "dem", "hear"]]
        + [f"results/rg_{p}.json" for p in
           ["pc_broad_dry", "pc_dry_wet", "mac_ckd", "mac_dem", "mac_glom",
            "mac_hear", "mac_htn", "mac_str", "mac_t2d"]]
        + ["results/traj_main_table.tsv", "results/traj_all_pairs.tsv",
           "results/ms_numbers.json", "docs/DATA_PROVENANCE.md"]
    )
    missing = [r for r in required if not os.path.exists(os.path.join(ROOT, r))]
    if missing:
        for r in missing:
            print("MISSING:", r)
        sys.exit(f"FAILED: {len(missing)} required file(s) missing")
    print(f"PASSED: {len(listed)} files verified, required tree complete")

if __name__ == "__main__":
    main()

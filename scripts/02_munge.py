#!/usr/bin/env python3
"""02_munge.py - filter raw FinnGen R11 sumstats to the 1000G EUR QC SNP set.
Logic identical to the sandbox pipeline: keep rsids in ldsnps.txt with
af_alt > 0.01 and sebeta > 0, Z = round(beta/sebeta, 4), keep-first dedupe by SNP.
Output per trait: derived/{trait}.sumstats.tsv.gz  (columns: SNP A1 A2 Z N P FRQ)
Idempotent via derived/_manifest.tsv (trait, kept_rows, output_bytes)."""
import os, sys, gzip, glob, subprocess, gc, time
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

CHUNK = 300_000

def ensure_snps():
    """Decompress resources/ldsnps.txt.gz -> derived/ldsnps.txt once; return path."""
    txt = common.p("derived", "ldsnps.txt")
    if not os.path.exists(txt) or os.path.getsize(txt) == 0:
        common.log("decompressing SNP list ...")
        with gzip_open(common.p("resources", "ldsnps.txt.gz")) as r, open(txt, "w") as w:
            for ln in r:
                w.write(ln)
    return txt

def gzip_open(path, mode="rt"):
    return gzip.open(path, mode)

def load_snp_set(txt):
    with open(txt) as f:
        return frozenset(ln.strip() for ln in f)

def munge_one(raw_gz, out_gz, snps_txt, n_total, workdir, log=common.log):
    """One trait: raw FinnGen gz -> filtered deduped gz. Returns kept row count.
    Raises on any inconsistency; safe to re-run (overwrites outputs)."""
    snps = load_snp_set(snps_txt)
    log(f"  snp set: {len(snps):,}")
    t0 = time.time()
    out_txt = os.path.join(workdir, "munge_tmp.tsv")
    kept = total = 0
    with open(out_txt, "w") as w:
        w.write("SNP\tA1\tA2\tZ\tN\tP\tFRQ\n")
        for ch in pd.read_csv(raw_gz, sep="\t", chunksize=CHUNK,
                              usecols=["rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"]):
            total += len(ch)
            # tolerate malformed source rows (empty/tab-collapsed fields force
            # object dtype; coerce to numeric so dropna can remove them)
            for col in ["pval", "beta", "sebeta", "af_alt"]:
                if ch[col].dtype == object:
                    ch[col] = pd.to_numeric(ch[col], errors="coerce")
            ch = ch.dropna(subset=["pval", "beta", "sebeta", "af_alt"])
            ch = ch[(ch.sebeta > 0) & (ch.af_alt > 0.01)]
            m = ch[ch.rsids.isin(snps)]
            if len(m):
                z = (m.beta / m.sebeta).round(4)
                blk = pd.DataFrame({"SNP": m.rsids, "A1": m.alt, "A2": m.ref,
                                    "Z": z, "N": n_total, "P": m.pval, "FRQ": m.af_alt})
                w.write(blk.to_csv(sep="\t", index=False, header=False))
                kept += len(m)
            del ch
    del snps
    gc.collect()
    out_dd = os.path.join(workdir, "munge_tmp_dd.tsv")
    subprocess.run(["bash", "-c",
        f"awk -F'\\t' 'NR==1{{print;next}} !seen[$1]++' '{out_txt}' > '{out_dd}'"], check=True)
    nrows = int(subprocess.run(["bash", "-c", f"wc -l < '{out_dd}'"],
                               capture_output=True, text=True).stdout) - 1
    subprocess.run(["bash", "-c", f"gzip -c '{out_dd}' > '{out_gz}'"], check=True)
    os.remove(out_txt); os.remove(out_dd)
    log(f"  rows_in={total:,} kept={nrows:,} ({time.time()-t0:.0f}s)")
    return nrows

def load_manifest():
    d = {}
    mf = common.p("derived", "_manifest.tsv")
    if os.path.exists(mf):
        for ln in open(mf):
            t = ln.rstrip("\n").split("\t")
            if len(t) == 3 and t[0] != "trait":
                d[t[0]] = (int(t[1]), int(t[2]))
    return d

def save_manifest(man):
    with open(common.p("derived", "_manifest.tsv"), "w") as f:
        f.write("trait\tkept\tbytes\n")
        for t, (k, b) in man.items():
            f.write(f"{t}\t{k}\t{b}\n")

def main():
    eps = common.load_endpoints()
    txt = ensure_snps()
    man = load_manifest()
    workdir = common.p("derived")
    for ep in eps:
        trait, endpoint, N = ep["trait"], ep["endpoint"], ep["N_total"]
        raw = common.p("data", "raw", f"finngen_R11_{endpoint}.gz")
        if not (os.path.exists(raw) and os.path.getsize(raw) == ep["remote_bytes"]):
            sys.exit(f"FATAL: {raw} missing or wrong size. Run with download first: bash run_all.sh")
        out = common.p("derived", f"{trait}.sumstats.tsv.gz")
        if trait in man and os.path.exists(out) and os.path.getsize(out) == man[trait][1]:
            common.log(trait, "verified done, skip")
            continue
        if os.path.exists(out):
            os.remove(out)
        common.log(trait, f"({endpoint}) munging ...")
        kept = munge_one(raw, out, txt, ep["N_total"], workdir)
        if kept < 5_000_000:
            sys.exit(f"FATAL: {trait} kept only {kept:,} SNPs (expected ~7-8M). "
                     "Do not continue; send back logs/ and results/metadata/.")
        man[trait] = (kept, os.path.getsize(out))
        save_manifest(man)
    common.log("MUNGE-ALL-DONE")

if __name__ == "__main__":
    main()

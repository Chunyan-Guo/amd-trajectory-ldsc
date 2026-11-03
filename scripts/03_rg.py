#!/usr/bin/env python3
"""03_rg.py - cross-trait LD score regression (estimator identical to the validated
sandbox implementation: per-chromosome weighted LS on (1, L2) with w=1/max(L2,1),
rg = slope_xy / sqrt(slope_a * slope_b), LOCO jackknife over 22 chromosomes).

Phases:
  main  - per-trait QC on V1 frame (qc_{trait}.json) + V1 (MAF>=0.05) and V3
          (V1 minus MHC chr6:25-34Mb) for all 9 pairs (rg_{pair}.json)
  v2    - V2 (MAF>=0.01) for all 9 pairs
Idempotent: completed (pair, variant) results are cached in their JSON files.
"""
import os, sys, glob, re, json, gc
import pandas as pd
import numpy as np
from scipy import stats as sstats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common

LD_DIR = common.p("resources", "ldscores")
MUNGE_DIR = common.p("derived")
RES_DIR = common.p("results")
LDF = sorted(glob.glob(os.path.join(LD_DIR, "ld*.l2.ldscore.gz")),
             key=lambda pth: int(re.search(r"ld(\d+)\.l2", pth).group(1)))
MHC = (6, 25_000_000, 34_000_000)

PAIRS = [("pc_broad_dry", "amd", "dry"), ("pc_dry_wet", "dry", "wet")] + \
        [("mac_" + b, "mac", b) for b in ["ckd", "dem", "glom", "hear", "htn", "str", "t2d"]]
TRAITS = ["mac", "amd", "dry", "wet", "htn", "str", "t2d", "ckd", "glom", "dem", "hear"]

def load_trait(name):
    fp = common.p("derived", f"{name}.sumstats.tsv.gz")
    if not os.path.exists(fp):
        sys.exit(f"FATAL: {fp} missing - run stage 02 first")
    df = pd.read_csv(fp, sep="\t", usecols=["SNP", "Z", "N", "FRQ"],
                     dtype={"Z": np.float32, "N": np.float64, "FRQ": np.float32})
    return df.set_index("SNP")

CACHE = {}; ORDER = []
def trait(name):
    if name in CACHE:
        return CACHE[name]
    df = load_trait(name)
    CACHE[name] = df; ORDER.append(name)
    while len(CACHE) > 2:
        old = ORDER.pop(0)
        if old in CACHE:
            del CACHE[old]
    gc.collect()
    return CACHE[name]

def accum_keys():
    return dict(sw=0., swx=0., swx2=0., sa=0., sxa=0., sb=0., sxb=0., sz=0., sxz=0., n=0)

def maf_mask005(j):
    return (np.minimum(j.Fa, 1 - j.Fa) >= 0.05) & (np.minimum(j.Fb, 1 - j.Fb) >= 0.05)

def maf_mask001(j):
    return (np.minimum(j.Fa, 1 - j.Fa) >= 0.01) & (np.minimum(j.Fb, 1 - j.Fb) >= 0.01)

def nomhc(j):
    return ~((j.CHR == MHC[0]) & (j.BP >= MHC[1]) & (j.BP <= MHC[2]))

def pair_sums_multi(dfa, dfb, variants):
    za = dfa.rename(columns={"Z": "Za", "N": "Na", "FRQ": "Fa"})
    zb = dfb.rename(columns={"Z": "Zb", "N": "Nb", "FRQ": "Fb"})
    MS = {"V0": lambda j: np.ones(len(j), bool),
          "V1": maf_mask005,
          "V2": maf_mask001,
          "V3": lambda j: maf_mask005(j) & nomhc(j)}
    S = {v: {c: accum_keys() for c in range(1, 23)} for v in variants}
    for f in LDF:
        c = int(re.search(r"ld(\d+)\.l2", f).group(1))
        dd = pd.read_csv(f, sep="\t")
        j = dd.join(za, on="SNP", how="inner").join(zb, on="SNP", how="inner")
        del dd
        if len(j) == 0:
            continue
        l2 = j.L2.values.astype(float)
        pa = (j.Za.values ** 2); pb = (j.Zb.values ** 2); pz = (j.Za.values * j.Zb.values)
        w = 1 / np.clip(l2, 1, None)
        for v in variants:
            keep = MS[v](j)
            if not keep.any():
                continue
            s = S[v][c]
            s["sw"] += w[keep].sum(); s["swx"] += (w * l2)[keep].sum(); s["swx2"] += (w * l2 * l2)[keep].sum()
            s["sa"] += (w * pa)[keep].sum(); s["sxa"] += (w * l2 * pa)[keep].sum()
            s["sb"] += (w * pb)[keep].sum(); s["sxb"] += (w * l2 * pb)[keep].sum()
            s["sz"] += (w * pz)[keep].sum(); s["sxz"] += (w * l2 * pz)[keep].sum()
            s["n"] += int(keep.sum())
        del j
    return S

def agg(S, key, skip=None):
    return sum(v[key] for k, v in S.items() if k != skip)

def fit2(S, sy, sxy, skip=None):
    A = np.array([[agg(S, "sw", skip), agg(S, "swx", skip)],
                  [agg(S, "swx", skip), agg(S, "swx2", skip)]])
    b = np.array([agg(S, sy, skip), agg(S, sxy, skip)])
    return np.linalg.solve(A, b)

def jack(vals):
    vals = np.asarray(vals)
    return float(np.sqrt(21 / 22 * np.sum((vals - vals.mean()) ** 2)))

def summarize(S):
    b11 = fit2(S, "sa", "sxa"); b22 = fit2(S, "sb", "sxb"); bxy = fit2(S, "sz", "sxz")
    def rg(skip=None):
        a = fit2(S, "sa", "sxa", skip)[1]; b = fit2(S, "sb", "sxb", skip)[1]; c = fit2(S, "sz", "sxz", skip)[1]
        return c / np.sqrt(a * b)
    rgs = [rg(skip=c) for c in range(1, 23)]
    h1 = [fit2(S, "sa", "sxa", skip=c)[1] for c in range(1, 23)]
    h2 = [fit2(S, "sb", "sxb", skip=c)[1] for c in range(1, 23)]
    rg_full = rg()
    se = jack(rgs)
    z = rg_full / se if se > 0 else float("nan")
    return dict(rg=float(rg_full), se=se, z=float(z),
                p=float(2 * sstats.norm.sf(abs(z))),
                h2a=float(b11[1]), h2a_se=jack(h1),
                h2b=float(b22[1]), h2b_se=jack(h2),
                int_xy=float(bxy[0]), int_a=float(b11[0]), int_b=float(b22[0]),
                nsnp=int(sum(v["n"] for v in S.values())))

def single_qc(df):
    S = {c: accum_keys() for c in range(1, 23)}
    for f in LDF:
        c = int(re.search(r"ld(\d+)\.l2", f).group(1))
        dd = pd.read_csv(f, sep="\t")
        j = dd.join(df[["Z"]].rename(columns={"Z": "Za"}), on="SNP", how="inner")
        del dd
        if len(j) == 0:
            continue
        l2 = j.L2.values.astype(float); pa = (j.Za.values ** 2)
        w = 1 / np.clip(l2, 1, None)
        s = S[c]
        s["sw"] += w.sum(); s["swx"] += (w * l2).sum(); s["swx2"] += (w * l2 * l2).sum()
        s["sa"] += (w * pa).sum(); s["sxa"] += (w * l2 * pa).sum(); s["n"] += len(j)
        del j
    h1 = [fit2(S, "sa", "sxa", skip=c)[1] for c in range(1, 23)]
    b11 = fit2(S, "sa", "sxa")
    z2 = df.Z.values.astype(float) ** 2
    return dict(nsnp=int(sum(v["n"] for v in S.values())),
                nsnp_raw=int(len(df)),
                mean_chi2=float(z2.mean()), lambda_gc=float(np.median(z2) / 0.4549),
                h2=float(b11[1]), h2_se=jack(h1), int_h2=float(b11[0]))

def run_pair(label, a, b, variants):
    out = common.p("results", f"rg_{label}.json")
    res = {}
    if os.path.exists(out):
        res = json.load(open(out))
    todo = [v for v in variants if v not in res]
    if not todo:
        common.log(label, "cached")
        return res
    dfa = trait(a); dfb = trait(b)
    S = pair_sums_multi(dfa, dfb, todo)
    for v in todo:
        res[v] = summarize(S[v])
        common.log(label, v, "rg=%.4f se=%.4f p=%.3g nsnp=%d"
                   % (res[v]["rg"], res[v]["se"], res[v]["p"], res[v]["nsnp"]))
    common.write_json(res, "results", f"rg_{label}.json")
    return res

def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else "main"
    if len(LDF) != 22:
        sys.exit(f"FATAL: expected 22 LD files in {LD_DIR}, found {len(LDF)}")
    if phase == "main":
        for t in TRAITS:
            qf = common.p("results", f"qc_{t}.json")
            if os.path.exists(qf):
                continue
            df = trait(t)
            maf = np.minimum(df.FRQ, 1 - df.FRQ)
            res = single_qc(df[maf >= 0.05])
            common.write_json(res, "results", f"qc_{t}.json")
            common.log("QC", t, res)
        for label, a, b in PAIRS:
            run_pair(label, a, b, ["V1", "V3"])
    elif phase == "v2":
        for label, a, b in PAIRS:
            run_pair(label, a, b, ["V2"])
    else:
        sys.exit(f"unknown phase {phase}")
    common.log("PHASE", phase, "DONE")

if __name__ == "__main__":
    main()

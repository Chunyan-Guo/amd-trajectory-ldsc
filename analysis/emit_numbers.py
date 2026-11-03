#!/usr/bin/env python3
"""Assemble every number the manuscript needs into ms_numbers.json.
Sources (all recomputed, verified):
  rg_{label}.json  - per-pair V1/V2/V3 results (main phase)
  qc_{trait}.json  - per-trait QC on V1 frame
  traj_main_table.tsv / traj_all_pairs.tsv - trajectory layer
Fails loudly if any required file/key is missing: no fabricated fallbacks.
"""
import json, os, glob, sys
import pandas as pd

D = os.environ.get("AMD_RESULTS", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"))
N = {}

def need(path):
    if not os.path.exists(path):
        sys.exit(f"MISSING {path}")
    return json.load(open(path))

def f3(x):  # 3 decimals
    return f"{x:.3f}"

def f2(x):
    return f"{x:.2f}"

def pfmt(p):
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m} x 10^{int(e)}"
    return f"{p:.4f}".rstrip("0").rstrip(".")

SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
def phtml(p):
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m} × 10{str(int(e)).translate(SUP)}"
    return f"{p:.4f}".rstrip("0").rstrip(".")

def ci(rg, se):
    return f"{rg - 1.96*se:.2f} to {rg + 1.96*se:.2f}"

PAIRS7 = ["mac_hear", "mac_t2d", "mac_htn", "mac_str", "mac_dem", "mac_ckd", "mac_glom"]
ENDPOINT_LABEL = {
 "mac_hear": "Hearing loss (H8_OTHHEARINGLOSS)",
 "mac_t2d":  "Type 2 diabetes (T2D)",
 "mac_htn":  "Hypertension (I9_HYPTENS)",
 "mac_str":  "Stroke (C_STROKE)",
 "mac_dem":  "Dementia (F5_DEMENTIA)",
 "mac_ckd":  "Chronic kidney disease (N14_CHRONKIDNEYDIS)",
 "mac_glom": "Glomerular disease (N14_GLOMERULAR)",
}
ENDPOINT_CODE = {
 "mac_hear": "H8_OTHHEARINGLOSS", "mac_t2d": "T2D", "mac_htn": "I9_HYPTENS",
 "mac_str": "C_STROKE", "mac_dem": "F5_DEMENTIA", "mac_ckd": "N14_CHRONKIDNEYDIS",
 "mac_glom": "N14_GLOMERULAR",
}
CASES = {  # FinnGen R11 manifest, verified
 "mac": (18777, 413784), "amd": (11023, 419198), "dry": (7589, 298486), "wet": (5890, 300152),
 "htn": (137312, 316345), "str": (47841, 328079), "t2d": (71728, 369007),
 "ckd": (11265, 436208), "glom": (6764, 446969), "dem": (20667, 428053), "hear": (10961, 397865),
}

def verdict(minrg, bonf):
    if bonf:
        return "Positive (Bonferroni)"
    if minrg <= 0.30:
        return "Excludes rg >= 0.3"
    if minrg <= 0.40:
        return "Borderline"
    return "Cannot exclude"

# ---- genetic layer ----
ps7 = {}
for lab in PAIRS7:
    r = need(f"{D}/rg_{lab}.json")["V1"]
    se = r["se"]
    minrg = 2.80 * se
    bonf = r["p"] < 0.0071
    ps7[lab] = dict(rg=r["rg"], se=se, p=r["p"], nsnp=r["nsnp"], int_xy=r["int_xy"],
                    h2b=r["h2b"], h2b_se=r["h2b_se"], ci=ci(r["rg"], se),
                    rg3=f3(r["rg"]), se3=f3(se), minrg=f2(minrg), minrg_f=minrg,
                    pstr=phtml(r["p"]), bonf=bonf, verdict=verdict(minrg, bonf),
                    code=ENDPOINT_CODE[lab])
# BH q across the 7 tests
from statsmodels.stats.multitest import multipletests
pv = [ps7[l]["p"] for l in PAIRS7]
_, q7, _, _ = multipletests(pv, method="fdr_bh")
for l, q in zip(PAIRS7, q7):
    ps7[l]["q"] = float(q)
    ps7[l]["q3"] = f"{q:.3f}"
    ps7[l]["qstr"] = phtml(q) if q < 1e-3 else f"{q:.3f}"

N["pairs7"] = ps7
N["h2_mac"] = need(f"{D}/qc_mac.json")["h2"]
N["h2_mac_se"] = need(f"{D}/qc_mac.json")["h2_se"]
N["h2_mac_3"] = f3(N["h2_mac"])
N["h2_mac_se_3"] = f3(N["h2_mac_se"])

# positive controls
pc1 = need(f"{D}/rg_pc_broad_dry.json")["V1"]
pc2 = need(f"{D}/rg_pc_dry_wet.json")["V1"]
N["pc1"] = dict(rg3=f3(pc1["rg"]), se3=f3(pc1["se"]), pstr=phtml(pc1["p"]), rg=pc1["rg"], se=pc1["se"], p=pc1["p"], ci=ci(pc1["rg"], pc1["se"]))
N["pc2"] = dict(rg3=f3(pc2["rg"]), se3=f3(pc2["se"]), pstr=phtml(pc2["p"]), rg=pc2["rg"], se=pc2["se"], p=pc2["p"], ci=ci(pc2["rg"], pc2["se"]))

# h2 range across the 11 traits (for Results text)
h2vals = sorted(v["h2f"] for v in (dict(h2f=need(f"{D}/qc_{t}.json")["h2"]) for t in
              ["mac", "amd", "dry", "wet", "htn", "str", "t2d", "ckd", "glom", "dem", "hear"]))
N["qc_range"] = dict(h2_lo=f3(h2vals[0]), h2_hi=f3(h2vals[-1]))

# qc table (S3)
qc = {}
for t in ["mac", "amd", "dry", "wet", "htn", "str", "t2d", "ckd", "glom", "dem", "hear"]:
    q = need(f"{D}/qc_{t}.json")
    qc[t] = dict(nsnp=q["nsnp"], nsnp_raw=q["nsnp_raw"], mean_chi2=round(q["mean_chi2"], 4),
                 lambda_gc=round(q["lambda_gc"], 4), h2=f3(q["h2"]), h2_se=f3(q["h2_se"]),
                 h2f=q["h2"], int_h2=f3(q["int_h2"]))
    qc[t]["cases"], qc[t]["controls"] = CASES[t]
N["qc"] = qc

# sensitivity V2/V3
sens = {}
for lab in PAIRS7:
    r = need(f"{D}/rg_{lab}.json")
    sens[lab] = {}
    for v in ("V1", "V2", "V3"):
        if v not in r:
            sys.exit(f"MISSING {v} in rg_{lab}.json")
        sens[lab][v] = dict(rg3=f3(r[v]["rg"]), se3=f3(r[v]["se"]), pstr=phtml(r[v]["p"]))
N["sens"] = sens
for lab in ("pc_broad_dry", "pc_dry_wet"):
    r = need(f"{D}/rg_{lab}.json")
    sens[lab] = {v: dict(rg3=f3(r[v]["rg"]), se3=f3(r[v]["se"]), pstr=phtml(r[v]["p"]))
                 for v in ("V1", "V2", "V3") if v in r}

# ---- trajectory layer ----
main = pd.read_csv(f"{D}/traj_main_table.tsv", sep="\t")
allp = pd.read_csv(f"{D}/traj_all_pairs.tsv", sep="\t")
sig = allp[allp["p_dir_fmt"].notna()]
N["traj"] = dict(n_records=int(len(allp)), n_pairs=int(allp["pair_id"].nunique()),
                 n_sig=int(len(sig)),
                 n_before=int((sig["relative_direction"] == "disease before H35").sum()),
                 n_after=int((sig["relative_direction"] == "H35 before disease").sum()),
                 n_rr1=int(len(main)),
                 n_rr1_nonobs=int((~main["is_obstetric_perinatal"].astype(bool)).sum()),
                 n_rr_lt1=int((sig["RR"] < 1).sum()))

def trow(r, direction):
    nm = str(r["name"])
    return dict(direction=direction, code=r["disease"], name=nm, n=int(r["counts"]),
                rr=f2(r["RR"]), pstr=r["p_dir_fmt"],
                q=(f"{r['q_fdr']:.1e}".replace("e", " x 10^") if pd.notna(r["q_fdr"]) else ""),
                age=f"{r['age_years']:.1f}", pf=(f"{r['pct_female']:.0f}" if pd.notna(r["pct_female"]) else "NA"),
                days=f"{r['CODE_DIFF_DAYS']:.0f}",
                obstetric=bool(r["is_obstetric_perinatal"]))
rows = []
for _, r in main[main.relative_direction == "disease before H35"].sort_values("RR", ascending=False).iterrows():
    rows.append(trow(r, "before"))
for _, r in main[main.relative_direction == "H35 before disease"].sort_values("RR", ascending=False).iterrows():
    rows.append(trow(r, "after"))
N["table1_rows"] = rows

json.dump(N, open(f"{D}/ms_numbers.json", "w"), indent=1)
print("ms_numbers.json written")
print("keys:", sorted(N.keys()))
print("hear:", ps7["mac_hear"]["rg3"], ps7["mac_hear"]["se3"], ps7["mac_hear"]["pstr"], ps7["mac_hear"]["ci"])
print("pc1:", N["pc1"]["rg3"], "+-", N["pc1"]["se3"], " pc2:", N["pc2"]["rg3"], "+-", N["pc2"]["se3"])

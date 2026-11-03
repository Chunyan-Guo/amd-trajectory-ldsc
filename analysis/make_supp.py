#!/usr/bin/env python3
"""Build Supplementary Tables S1-S4 (one xlsx, 4 sheets).
S1: all 421 ordered H35 records from the DTB aggregate export + exclusion/flag reasons
S2: age/sex audit of the 199 directioned records
S3: per-trait and per-pair LDSC QC
S4: sensitivity analyses (MAF>=0.01; minus MHC)
Fails loudly on any missing input.
"""
import json, os, sys
import pandas as pd

D = os.environ.get("AMD_RESULTS", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"))
OUT = os.path.join(D, "supplement_S1_S4.xlsx")

def need(p):
    if not os.path.exists(p):
        sys.exit(f"MISSING {p}")
    return p

# ---------------- S1 ----------------
allp = pd.read_csv(need(f"{D}/traj_all_pairs.tsv"), sep="\t")
assert len(allp) == 421, f"S1 expects 421 records, got {len(allp)}"

def p2f(s):
    s = str(s)
    if s.startswith("<"): return float(s[1:])
    try: return float(s)
    except Exception: return float("nan")

allp["p_dir_numeric"] = allp["p_dir_fmt"].map(p2f)
allp["binom_floor"] = 2 * 0.5 ** allp["counts"].fillna(10**9)
subfloor = allp["p_dir_numeric"] < allp["binom_floor"]
allp["published_P_below_exact_binomial_floor"] = ["yes" if b else "" for b in subfloor]
s1 = allp.rename(columns={
    "pair_id": "Pair ID", "disease": "ICD-10 code", "name": "Diagnosis name",
    "relative_direction": "Direction relative to H35", "RR": "Relative risk (5-year co-occurrence)",
    "counts": "Patients, n", "pct_female": "Female, %", "age_years": "Mean age at first diagnosis, y",
    "CODE_DIFF_DAYS": "Median days between diagnoses", "death5y_n": "Deaths within 5 y, n",
    "p_5yr_fmt": "5-year co-occurrence P", "p_dir_fmt": "Directionality P (as published)",
    "q_fdr": "FDR q (Benjamini-Hochberg, this study)", "flag": "Flag / exclusion reason"})
s1 = s1[["Pair ID", "ICD-10 code", "Diagnosis name", "Direction relative to H35",
         "Relative risk (5-year co-occurrence)", "Patients, n", "Female, %",
         "Mean age at first diagnosis, y", "Median days between diagnoses", "Deaths within 5 y, n",
         "5-year co-occurrence P", "Directionality P (as published)",
         "FDR q (Benjamini-Hochberg, this study)", "published_P_below_exact_binomial_floor",
         "Flag / exclusion reason"]]
n_subfloor = int(subfloor.sum())
print(f"S1: {len(s1)} records; {n_subfloor} flagged sub-floor P")

# ---------------- S2 ----------------
sig = allp[allp["q_fdr"].notna()].copy()
assert len(sig) == 199, f"S2 expects 199 directioned records, got {len(sig)}"
s2 = sig.rename(columns={
    "disease": "ICD-10 code", "name": "Diagnosis name",
    "relative_direction": "Direction relative to H35", "RR": "Relative risk",
    "counts": "Patients, n", "age_years": "Mean age at first diagnosis, y",
    "pct_female": "Female, %", "flag": "Audit flag"})
sex_specific = s2["ICD-10 code"].str[0].isin(["N", "O"]) | s2["ICD-10 code"].isin(
    ["C62", "D06", "D24", "D25", "C50", "D27", "D28", "N60", "N87", "N92"])
s2["Sex- or age-structure related code"] = ["yes" if x else "" for x in sex_specific]
s2 = s2[["ICD-10 code", "Diagnosis name", "Direction relative to H35", "Relative risk",
         "Patients, n", "Mean age at first diagnosis, y", "Female, %",
         "Sex- or age-structure related code", "Audit flag"]]
print(f"S2: {len(s2)} directioned records audited")

# ---------------- S3 / S4 ----------------
NUM = json.load(open(need(f"{D}/ms_numbers.json")))
TRAITS = ["mac", "amd", "dry", "wet", "htn", "str", "t2d", "ckd", "glom", "dem", "hear"]
TRAIT_CODE = {"mac": "H7_MACULADEGEN", "amd": "H7_AMD", "dry": "DRY_AMD", "wet": "WET_AMD",
              "htn": "I9_HYPTENS", "str": "C_STROKE", "t2d": "T2D",
              "ckd": "N14_CHRONKIDNEYDIS", "glom": "N14_GLOMERULAR",
              "dem": "F5_DEMENTIA", "hear": "H8_OTHHEARINGLOSS"}
rows = []
for t in TRAITS:
    q = NUM["qc"][t]
    rows.append({"Trait": t, "FinnGen R11 endpoint": TRAIT_CODE[t],
                 "Cases": q["cases"], "Controls": q["controls"],
                 "SNPs after merge+MAF>=0.05": q["nsnp"], "SNPs before QC": q["nsnp_raw"],
                 "Mean chi2": q["mean_chi2"], "lambda GC": q["lambda_gc"],
                 "h2 (observed)": q["h2"], "h2 SE": q["h2_se"], "h2 intercept": q["int_h2"]})
s3a = pd.DataFrame(rows)

pairs_rows = []
PAIR_LABEL = {"pc_broad_dry": "H7_AMD x DRY_AMD (positive control 1)",
              "pc_dry_wet": "DRY_AMD x WET_AMD (positive control 2)",
              "mac_htn": "H7_MACULADEGEN x I9_HYPTENS", "mac_str": "H7_MACULADEGEN x C_STROKE",
              "mac_t2d": "H7_MACULADEGEN x T2D", "mac_ckd": "H7_MACULADEGEN x N14_CHRONKIDNEYDIS",
              "mac_glom": "H7_MACULADEGEN x N14_GLOMERULAR", "mac_dem": "H7_MACULADEGEN x F5_DEMENTIA",
              "mac_hear": "H7_MACULADEGEN x H8_OTHHEARINGLOSS"}
for lab in PAIR_LABEL:
    r = json.load(open(need(f"{D}/rg_{lab}.json")))["V1"]
    pairs_rows.append({"Pair": PAIR_LABEL[lab], "rg": round(r["rg"], 4), "SE": round(r["se"], 4),
                       "P": r["p"], "n SNP": r["nsnp"], "gcov intercept": round(r["int_xy"], 4),
                       "h2 trait A": round(r["h2a"], 4), "h2 trait B": round(r["h2b"], 4)})
s3b = pd.DataFrame(pairs_rows)

sens_rows = []
for lab in PAIR_LABEL:
    r = json.load(open(need(f"{D}/rg_{lab}.json")))
    for v, desc in [("V1", "Primary (MAF >= 0.05)"), ("V2", "MAF >= 0.01"), ("V3", "MAF >= 0.05, minus MHC")]:
        if v not in r: sys.exit(f"MISSING {v} in rg_{lab}.json")
        sens_rows.append({"Pair": PAIR_LABEL[lab], "Analysis": desc,
                          "rg": round(r[v]["rg"], 4), "SE": round(r[v]["se"], 4), "P": r[v]["p"]})
s4 = pd.DataFrame(sens_rows)

with pd.ExcelWriter(OUT, engine="openpyxl") as xw:
    s1.to_excel(xw, sheet_name="S1_all_421_records", index=False)
    s2.to_excel(xw, sheet_name="S2_age_sex_audit", index=False)
    s3a.to_excel(xw, sheet_name="S3_ldsc_qc_traits", index=False)
    s3b.to_excel(xw, sheet_name="S3_ldsc_qc_pairs", index=False)
    s4.to_excel(xw, sheet_name="S4_sensitivity", index=False)
print("WROTE", OUT)

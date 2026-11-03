#!/usr/bin/env python3
"""Rebuild the trajectory-layer tables from the shipped DTB aggregate export.

Reads:  data/dtb_supp1.tsv  (Danish Disease Trajectory Browser aggregate export,
        Siggaard et al., Nature Communications 2020; shipped with this package)
Writes: results/traj_all_pairs.tsv   (all 421 directed records involving H35)
        results/traj_main_table.tsv  (Table 1 content: directioned, RR>1, q<0.05)

Filters used, exactly as in the manuscript: records involving ICD-10 H35;
directioned subset = records with direction_yes_no == 1; Benjamini-Hochberg FDR
across the directioned subset; main table keeps RR > 1 and q < 0.05.
The DTB browser itself exposes the same records (diagnosis pair query on H35);
this script reproduces the derived tables deterministically from the export.
"""
import os, sys
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from icd10_dict import ICD10

dtb = pd.read_csv(os.path.join(ROOT, "data", "dtb_supp1.tsv"), sep="\t")
h35 = dtb[(dtb.D1 == "H35") | (dtb.D2 == "H35")].copy()

h35["disease"] = np.where(h35.D1 == "H35", h35.D2, h35.D1)
h35["relative_direction"] = np.where(h35.D1 == "H35", "H35 before disease", "disease before H35")
h35["name"] = h35.disease.map(ICD10)
h35["chapter"] = h35.disease.str[0]
h35["is_obstetric_perinatal"] = h35.chapter.isin(["O", "P"])
h35["age_years"] = (h35.AGE_AT_DISEASE / 365.25).round(1)
h35["pct_female"] = (100 * h35.female_counts / h35.counts).round(1)
h35["death5y_n"] = h35.death_counts

pairs = h35.apply(lambda r: tuple(sorted([r.D1, r.D2])), axis=1)
h35["pair_id"] = pairs.map(lambda t: t[0] + "_" + t[1])

sig = h35[h35.direction_yes_no == 1].copy()
rej, q, _, _ = multipletests(sig["p.value.direction"].clip(lower=1e-300), method="fdr_bh")
sig["q_fdr"] = q
h35 = h35.merge(sig[["pair_id", "relative_direction", "q_fdr"]],
                on=["pair_id", "relative_direction"], how="left")

def fmt_p(p):
    if pd.isna(p): return ""
    if p == 0: return "<1e-300"
    return "%.2e" % p

h35["p_dir_fmt"] = h35["p.value.direction"].map(fmt_p)
h35["p_5yr_fmt"] = h35.p5years.map(fmt_p)

def flag(r):
    fs = []
    if r.direction_yes_no != 1: fs.append("no significant directionality (reverse record)")
    if r.direction_yes_no == 1 and r.RR <= 1: fs.append("RR<1 (negative association; competing risk / surveillance / age-structure artefact)")
    if r.is_obstetric_perinatal: fs.append("obstetric/perinatal code; age-structure contamination (mean age %.1f y)" % (r.age_years if pd.notna(r.age_years) else float("nan")))
    return "; ".join(fs)
h35["flag"] = h35.apply(flag, axis=1)

res = os.path.join(ROOT, "results")
cols = ["pair_id", "disease", "name", "relative_direction", "RR", "counts", "pct_female",
        "age_years", "CODE_DIFF_DAYS", "death5y_n", "p_5yr_fmt", "p_dir_fmt", "q_fdr", "flag"]
h35.sort_values(["pair_id", "relative_direction"])[cols].to_csv(os.path.join(res, "traj_all_pairs.tsv"), sep="\t", index=False)

main = sig[(sig.RR > 1) & (sig.q_fdr < 0.05)].copy()
main["name"] = main.disease.map(ICD10)
main["p_dir_fmt"] = main["p.value.direction"].map(fmt_p)
main = main.sort_values(["relative_direction", "RR"], ascending=[True, False])
mcols = ["disease", "name", "relative_direction", "RR", "counts", "pct_female", "age_years",
         "CODE_DIFF_DAYS", "p_dir_fmt", "q_fdr", "is_obstetric_perinatal"]
main[mcols].to_csv(os.path.join(res, "traj_main_table.tsv"), sep="\t", index=False)
print("traj_all_pairs.tsv:", len(h35), "records; traj_main_table.tsv:", len(main), "records")

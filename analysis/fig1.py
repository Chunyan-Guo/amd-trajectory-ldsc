#!/usr/bin/env python3
"""Figure 1: disease trajectories to/from H35 (retinal degeneration), RR>1, FDR q<0.05.
Strict inch-grid layout, fixed-offset panel letters, log x-axis, point area ~ n.
"""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator, NullFormatter

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"

import os
D = os.environ.get("AMD_RESULTS", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"))
T = pd.read_csv(os.path.join(D, "traj_main_table.tsv"), sep="\t")

SHORT = {
 "H54": "Visual impairment including blindness",
 "H57": "Other disorders of eye and adnexa",
 "N08": "Glomerular disorders in other diseases",
 "E16": "Other pancreatic endocrine disorders",
 "M32": "Systemic lupus erythematosus",
 "E10": "Type 1 diabetes mellitus",
 "H10": "Conjunctivitis",
 "O24": "Diabetes mellitus in pregnancy",
 "M31": "Other necrotizing vasculopathies",
 "I65": "Occlusion/stenosis of precerebral arteries",
 "E11": "Type 2 diabetes mellitus",
 "H81": "Disorders of vestibular function",
 "M35": "Other systemic connective tissue involvement",
 "I10": "Essential (primary) hypertension",
 "I20": "Angina pectoris",
 "G45": "Transient cerebral ischaemic attacks",
 "E78": "Disorders of lipoprotein metabolism",
}
T["short"] = T.disease.map(SHORT)
assert T.short.notna().all()

before = T[T.relative_direction == "disease before H35"].sort_values("RR", ascending=False).reset_index(drop=True)
after  = T[T.relative_direction == "H35 before disease"].sort_values("RR", ascending=False).reset_index(drop=True)

# ---- strict inch grid ----
W = 7.2
LBL_W = 2.85
PLOT_W = 3.20
TOP_M, BOT_M = 0.40, 0.60
ROW = 0.315
GAP = 0.85          # legend lives in the gap

hA = len(before) * ROW
hB = len(after) * ROW
H = TOP_M + hA + GAP + hB + BOT_M

fig = plt.figure(figsize=(W, H))
axw = PLOT_W / W
axA = fig.add_axes([LBL_W / W, (H - TOP_M - hA) / H, axw, hA / H])
axB = fig.add_axes([LBL_W / W, BOT_M / H, axw, hB / H])

XLO, XHI = 1.0, 6.0
TICKS = [1, 1.5, 2, 3, 4]
C_A = "#0072B2"
C_B = "#D55E00"
C_GREY = "#999999"

def nsize(n):
    return 22 + 400 * n / T.counts.max()

def draw(ax, df, color, title, show_xlabels):
    ax.set_xscale("log")
    ax.set_xlim(XLO, XHI)
    ax.set_ylim(-1.0, len(df) - 0.4)
    ax.invert_yaxis()
    ax.set_yticks(range(len(df)))
    labs = []
    for r in df.itertuples():
        s = f"{r.disease}  {r.short}"
        if r.disease == "O24":
            s += " \u2020"
        labs.append(s)
    ax.set_yticklabels(labs, fontsize=6.9)
    ax.set_xticks(TICKS)
    if show_xlabels:
        ax.set_xticklabels([str(t) for t in TICKS], fontsize=7.0)
    else:
        ax.set_xticklabels([])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis="both", length=2.5, width=0.6, labelsize=7.0)
    ax.grid(axis="x", which="major", color="#DDDDDD", lw=0.5, zorder=0)
    ax.axvline(1.0, color="#444444", lw=0.8, zorder=1)
    ax.text(5.85, -0.68, "RR (n)", ha="right", va="center", fontsize=6.9,
            fontweight="bold", color="#222222")
    for i, r in enumerate(df.itertuples()):
        c = C_GREY if r.disease == "O24" else color
        ax.scatter([r.RR], [i], s=nsize(r.counts), c=c, alpha=0.85,
                   edgecolors="white", linewidths=0.5, zorder=3)
        ax.text(5.85, i, f"{r.RR:.2f}  (n={int(r.counts):,})",
                ha="right", va="center", fontsize=6.6, color="#222222")
    ax.set_title(title, loc="left", fontsize=8.4, fontweight="bold", pad=6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(0.7)

draw(axA, before, C_A, "Diagnoses preceding retinal degeneration (H35)", show_xlabels=False)
draw(axB, after,  C_B, "Diagnoses following retinal degeneration (H35)", show_xlabels=True)
axB.set_xlabel("Relative risk (log scale)", fontsize=7.6, labelpad=3)

for ax, letter in ((axA, "A"), (axB, "B")):
    pos = ax.get_position()
    fig.text(pos.x0 - 0.20 / W, pos.y1 + 0.10 / H, letter,
             fontsize=10.5, fontweight="bold", va="bottom", ha="left")

# horizontal size legend centered in the inter-panel gap
leg_n = [100, 1000, 10000]
handles = [Line2D([], [], marker="o", ls="", markerfacecolor="#BBBBBB",
                  markeredgecolor="white", markersize=np.sqrt(nsize(n)),
                  label=f"{n:,}") for n in leg_n]
gap_c = (BOT_M + hB + GAP * 0.50) / H
fig.legend(handles=handles, ncol=3, columnspacing=1.4, handletextpad=0.55,
           loc="center", bbox_to_anchor=(0.545, gap_c),
           fontsize=6.6, frameon=False, borderpad=0.0, labelspacing=0.0)
fig.text(0.545, gap_c + 0.32 / H, "Patients (n)", fontsize=6.8,
         fontweight="bold", ha="center", va="center")

fig.text(0.015, 0.018,
         "\u2020 Obstetric code; interpreted as age-structure contamination (mean age 31.1 y, 100% female).\n"
         "95% CIs for RRs are not released in the aggregate export; point size is proportional to patient count.",
         fontsize=6.2, color="#555555", ha="left", va="bottom")

out1 = os.path.join(D, "fig1_trajectories.png")
out2 = os.path.join(D, "fig1_trajectories.pdf")
fig.savefig(out1, dpi=300, facecolor="white")
fig.savefig(out2, facecolor="white")
print("saved", W, round(H, 2))

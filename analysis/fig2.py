#!/usr/bin/env python3
"""Figure 2 - genetic-correlation forest plot (mac vs 7 systemic endpoints).
Strict inch-grid layout; fixed-offset elements; all-English; 300 dpi.
Reads ms_numbers.json (recompute pipeline outputs). No fabricated values.
"""
import json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = os.environ.get("AMD_RESULTS", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results"))
NUM = json.load(open(f"{D}/ms_numbers.json"))
P7 = NUM["pairs7"]

LBL = {"mac_hear": "Hearing loss", "mac_t2d": "Type 2 diabetes", "mac_htn": "Hypertension",
       "mac_str": "Stroke", "mac_dem": "Dementia", "mac_ckd": "Chronic kidney disease",
       "mac_glom": "Glomerular disease"}

order = sorted(P7.keys(), key=lambda k: P7[k]["rg"], reverse=True)
rows = [(LBL[k], P7[k]["code"], P7[k]["rg"], P7[k]["se"], P7[k]["bonf"], P7[k]["ci"], P7[k]["pstr"]) for k in order]

W, H = 7.2, 4.3
ML, LBL_W, PLOT_W, TXT_W = 0.12, 2.30, 2.55, 2.05   # inches
fig = plt.figure(figsize=(W, H))
AXY0, AXH = 0.21, 0.63
ax = fig.add_axes([ (ML+LBL_W)/W, AXY0, PLOT_W/W, AXH ])
ax.grid(False)

ys = np.arange(len(rows))[::-1]
C_POS, C_NULL = "#D55E00", "#0072B2"   # Okabe-Ito vermillion / blue
for y, (name, code, rg, se, bonf, ci, pstr) in zip(ys, rows):
    lo, hi = rg - 1.96*se, rg + 1.96*se
    c = C_POS if bonf else C_NULL
    ax.plot([lo, hi], [y, y], color=c, lw=1.4, zorder=2, solid_capstyle="round")
    ax.plot([lo, lo], [y-0.16, y+0.16], color=c, lw=1.2, zorder=2)
    ax.plot([hi, hi], [y-0.16, y+0.16], color=c, lw=1.2, zorder=2)
    ax.scatter([rg], [y], s=58 if bonf else 46, marker="D" if bonf else "o",
               facecolor=c if bonf else "white", edgecolor=c, linewidth=1.4, zorder=3)
ax.axvline(0, color="#555555", lw=1.0, ls=(0, (4, 3)), zorder=1)

allx = [r[2] for r in rows]
lo_min = min(r[2]-1.96*r[3] for r in rows); hi_max = max(r[2]+1.96*r[3] for r in rows)
x0 = np.floor((lo_min-0.05)*10)/10; x1 = np.ceil((hi_max+0.05)*10)/10
ax.set_xlim(x0, x1)
ax.set_ylim(-0.7, len(rows)-0.3)
step = 0.2
ax.set_xticks(np.arange(np.ceil(x0/step)*step, x1+1e-9, step))
ax.set_xticklabels([f"{v:.1f}" for v in np.arange(np.ceil(x0/step)*step, x1+1e-9, step)], fontsize=8)
ax.set_yticks([])
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="x", length=3, width=0.8)
ax.set_xlabel("Genetic correlation (rg) with macular degeneration", fontsize=9)

# left labels (name + FinnGen code, two lines)
axl = fig.add_axes([ML/W, AXY0, LBL_W/W, AXH]); axl.axis("off")
axl.set_xlim(0, 1); axl.set_ylim(-0.7, len(rows)-0.3)
for y, (name, code, *_x) in zip(ys, rows):
    axl.text(0.98, y+0.13, name, ha="right", va="center", fontsize=8.4)
    axl.text(0.98, y-0.20, code, ha="right", va="center", fontsize=6.6, color="#555555")

# right text column: rg (95% CI) and P
axt = fig.add_axes([(ML+LBL_W+PLOT_W+0.08)/W, AXY0, TXT_W/W, AXH]); axt.axis("off")
axt.set_xlim(0, 1); axt.set_ylim(-0.7, len(rows)-0.3)
axt.text(0, len(rows)-0.55, "rg (95% CI)", fontsize=8, fontweight="bold", va="center")
axt.text(0.66, len(rows)-0.55, "P", fontsize=8, fontweight="bold", va="center")
for y, (name, code, rg, se, bonf, ci, pstr) in zip(ys, rows):
    axt.text(0, y, ci, fontsize=7.6, va="center")
    axt.text(0.66, y, pstr, fontsize=7.6, va="center")

fig.text((ML+LBL_W-0.06)/W, 0.905, "Systemic endpoint (FinnGen R11)", fontsize=8, fontweight="bold", ha="right")

# footnote
fig.text(ML/W, 0.028,
         "Diamond: pair surviving Bonferroni correction (P < 0.0071). Whiskers: 95% CI (leave-one-chromosome-out jackknife).\n"
         "Dashed line: rg = 0. Full results including sensitivity analyses: Table 2 and Supplementary Tables S3-S4.",
         fontsize=7.2, color="#333333", va="bottom")

fig.savefig(f"{D}/fig2_rg_forest.png", dpi=300)
matplotlib.rcParams["pdf.fonttype"] = 42
fig.savefig(f"{D}/fig2_rg_forest.pdf")
print("WROTE fig2_rg_forest.png/.pdf", os.path.getsize(f"{D}/fig2_rg_forest.png"))

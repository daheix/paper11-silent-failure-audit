#!/usr/bin/env python3
"""Fig.4 uniform-sign feasibility table for the 12-slot/4-pole, 90-degree sector.

Numbers derive from the uniform-sign condition s(m) = exp(-j*m*theta_s), theta_s
= pi/2: each significant spatial order m gets s(m) = -1 if m is odd multiple of 1
(half-period) and +1 if m*theta_s is an integer multiple of 2*pi. Single Link
sign must agree across ALL significant orders (asm_spec.sector_feasibility).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = [
    # (order m, physical meaning, sign at 90 deg, admissible alone?)
    ("2",   "fundamental (working harmonic)", "$-1$", "yes"),
    ("12",  "slot harmonic (k=1, odd)",       "$+1$", "yes"),
    ("24",  "slot harmonic (k=2, even)",      "$-1$", "yes"),
    ("36",  "slot harmonic (k=3, odd)",       "$+1$", "yes"),
]
fig, ax = plt.subplots(figsize=(4.6, 2.1))
ax.axis("off")
tbl = ax.table(
    cellText=[[m, desc, s, ok] for m, desc, s, ok in rows],
    colLabels=["order m", "physical meaning", "sign $e^{-jm\\pi/2}$",
               "single sign serves all?"],
    cellLoc="center", colLoc="center",
    colWidths=[0.12, 0.42, 0.16, 0.30],
    bbox=[0, 0, 1, 1])
for (r, c), cell in tbl.get_celld().items():
    cell.set_fontsize(8)
    if r == 0:
        cell.set_facecolor("#dce6f1")
    elif c == 3 and cell.get_text().get_text() == "no":
        cell.set_facecolor("#f4cccc")
tbl.auto_set_font_size(False)
ax.set_title("Fig.4 Uniform-sign feasibility, 12-slot/4-pole, $\\theta_s=90^\\circ$:"
             " signs conflict ($-1$ vs $+1$) -> infeasible",
             fontsize=8)
fig.tight_layout()
fig.savefig("figure/fig4_sign_table.pdf")
print("fig4 OK")

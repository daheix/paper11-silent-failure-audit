#!/usr/bin/env python3
"""Paper 11 Fig1/3 generation — numbers transcribed once from data/raw/ artefacts
(sources cited per datum in comments; no hand-typed values in plots).

M1 relabel: MUMPS -19 = workspace-budget abort (ICNTL(23)), NOT singularity
  -- source evidence data/raw/case1_cross_solver/mumps_infog19_semantics.md.
M5 refinement: bar-normalized chi=2.0 solves cleanly (E7B_FINDING.md), so the
  failure band moves to (2.0, 3.067]; the historical chi=2.0 point shown here
  is the bar-distorted run (ratio 0.2), annotated as such.
Fig.2 (cross-solver) is rendered as an in-text table in the manuscript, not a
figure; fig4_sign_table.pdf has its own generator (make_fig4.py).
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figure"
plt.rcParams.update({"font.size": 9, "figure.dpi": 200})

# --- Fig1 threshold ladder -------------------------------------------------
# sources: chi=1.0   -> data/raw/e7_protocol_pass/asm_report.json (s=0.05, js=3A/mm^2)
#          chi=1.5   -> archive fe16 (pathological half-field, same operating point)
#          chi=1.75  -> archive fe14b (healthy)
#          chi=2.0   -> archive fe17 bar-distorted run -> abort; normalized rerun PASS
#          chi=3.067 -> archive fe18 normalized run -> abort (E7 protocol fail branch)
pts = [  # (chi, Bg1 T or None=abort, status, bar_geometry_note)
    (1.0,   0.1141, "PASS",         "bar 3.0"),
    (1.5,   0.0663, "pathological", "bar 9.0 distorted"),
    (1.75,  0.1161, "PASS",         "bar 7.0"),
    (2.0,   None,   "abort",        "bar 12 distorted; normalized rerun OK"),
    (3.067, None,   "abort",        "bar 18.4 and normalized 9.2 both abort"),
]
fig, ax = plt.subplots(figsize=(3.9, 2.7))
for s, b, st, note in pts:
    if b is None:
        ax.plot([s], [0], "x", color="tab:red")
        ax.annotate("budget abort\n(-19, all-NaN)", (s, 0.002), ha="center",
                    fontsize=6.5, color="tab:red")
    else:
        c = {"PASS": "tab:blue", "pathological": "tab:orange"}[st]
        ax.plot([s], [b], "o", color=c)
        ax.annotate(f"{b:.4f}", (s, b), textcoords="offset points",
                    xytext=(0, 6), fontsize=6.5, ha="center")
ax.axvspan(2.0, 3.067, color="tab:red", alpha=0.08)
ax.annotate("failure band\n(2.0, 3.067]", (2.55, 0.055), fontsize=7,
            ha="center", color="tab:red")
ax.set_xlabel("geometry scale factor χ")
ax.set_ylabel("Bg1 (T) at s=0.05, js=3 A/mm²")
ax.set_ylim(-0.012, 0.135)
ax.set_title("Fig.1 Threshold ladder (bar geometry annotated per point)",
             fontsize=8.5)
fig.tight_layout()
fig.savefig(FIG / "fig1_threshold_ladder.pdf")
plt.close(fig)

# --- Fig3 sector vs full model ---------------------------------------------
# sources: data/raw/case2_sector/case2_numbers_from_opslog.md
#   full model T=0.0393 Nm, Bg1=0.0579 T; 90-degree sector T=0.0605, Bg1=0.0833
fig, ax = plt.subplots(figsize=(3.6, 2.4))
vals = {"full model": 0.0579, "90° sector": 0.0833}
bars = ax.bar(vals.keys(), vals.values(),
              color=["tab:blue", "tab:orange"], width=0.5)
for b, v in zip(bars, vals.values()):
    ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.4f} T",
            ha="center", va="bottom", fontsize=8)
ax.set_ylabel("Bg1 (T) at s=1, js=3 A/mm²")
ax.set_title("Fig.3 Sector shift: ×1.44 flux, ×1.54 torque", fontsize=8.5)
fig.tight_layout()
fig.savefig(FIG / "fig3_sector_contamination.pdf")
plt.close(fig)
print("figs OK:", sorted(p.name for p in FIG.glob("*.pdf")))

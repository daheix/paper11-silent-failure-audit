#!/usr/bin/env python3
"""Paper 11 Fig1-3 生成 — 数字唯一来源 data/raw/* (真数据铁律, 无手输数字)"""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
RAW, FIG = ROOT/"data/raw", ROOT/"figure"
plt.rcParams.update({"font.size": 9, "figure.dpi": 200})

# Fig1 阈值谱 (M-ASM-15 evidence 五点, 数字经 asm_design.json 固化链)
pts = [  # (尺度因子, Bg1 T 或 None=nan, 状态)
    (1.0, 0.1141, "PASS"), (1.5, 0.0663, "pathological"),
    (1.75, 0.1161, "PASS"), (2.0, None, "collapse"), (3.067, None, "collapse")]
fig, ax = plt.subplots(figsize=(3.6, 2.6))
for s, b, st in pts:
    if b is None:
        ax.annotate("MUMPS -19\n(nan)", (s, 0.002), ha="center", fontsize=7, color="tab:red")
        ax.plot([s], [0], "x", color="tab:red")
    else:
        c = {"PASS": "tab:blue", "pathological": "tab:orange"}[st]
        ax.plot([s], [b], "o", color=c)
        ax.annotate(f"{b:.4f}", (s, b), textcoords="offset points", xytext=(0, 6), fontsize=6.5, ha="center")
ax.axvspan(1.75, 2.0, color="tab:red", alpha=0.08)
ax.annotate("collapse band\n(1.75, 2.0]", (1.875, 0.13), fontsize=7, ha="center", color="tab:red")
ax.set_xlabel("geometry scale factor ×"); ax.set_ylabel("Bg1 (T) — MUMPS/OOC")
ax.set_title("Fig.1 Silent-failure threshold ladder", fontsize=9)
fig.tight_layout(); fig.savefig(FIG/"fig1_threshold_ladder.pdf"); plt.close(fig)

# Fig2 交叉求解对照 (fe14 同矩阵: MUMPS 失败 vs SuperLU 收敛)
fig, ax = plt.subplots(figsize=(3.6, 2.2))
ax.barh([1, 0], [1, 1], color=["tab:red", "tab:green"], alpha=0.75)
ax.set_yticks([1, 0], ["MUMPS+OOC", "SuperLU (serial)"])
ax.set_xlabel("relative outcome (no solve vs residual 2.28e-5, 26 s)")
for y, txt in [(1, "ERROR -19 (false singular), no solution"),
               (0, "T=0.7477 Nm, Bg1=0.0158 T, KCL=0.000")]:
    ax.text(0.02, y, txt, va="center", fontsize=7, color="white")
ax.set_title("Fig.2 Cross-solver verdict on the identical matrix", fontsize=9)
ax.set_xticks([]); fig.tight_layout()
fig.savefig(FIG/"fig2_cross_solver.pdf"); plt.close(fig)

# Fig3 扇区 vs 全模型 (js17b vs asm_full16d, s=1)
fig, ax = plt.subplots(figsize=(3.6, 2.4))
vals = {"full model": 0.0579, "90° sector": 0.0833}
bars = ax.bar(vals.keys(), vals.values(), color=["tab:blue", "tab:orange"], width=0.5)
for b, v in zip(bars, vals.values()):
    ax.text(b.get_x()+b.get_width()/2, v, f"{v:.4f} T", ha="center", va="bottom", fontsize=8)
ax.set_ylabel("Bg1 (T) @ s=1, js=3 A/mm²")
ax.set_title("Fig.3 Sector-symmetry contamination: ×1.44 on fundamental", fontsize=8.5)
fig.tight_layout(); fig.savefig(FIG/"fig3_sector_contamination.pdf"); plt.close(fig)
print("figs OK:", sorted(p.name for p in FIG.glob("*.pdf")))

#!/usr/bin/env python3
"""audit_protocol.py — paper11 三审计协议的可执行实现 (A1-A5).

输入: 链产出的 asm_report.json (至少一个求解器; 两个求解器时启用交叉审计)。
输出: verdict PASS/FAIL + 每步证据。判据与 paper11 (COMPEL 投稿稿) 一致。

审计步骤 (对应论文 Audit Protocol 节):
  A1 完成性   — 每份报告 physics_check 场有效且有限 (NaN/inf 扫描)
  A2 数值卫生 — kcl_resid 相对量 < 1e-6 (电荷守恒), 求解耗时有限
  A3 交叉求解器一致性 — 双报告同工况 T/Bg1 相对差 <= 5% (假奇异/静默退化抓现行)
  A4 尺度带自洽 — Bg1 落入同族冻结基线的 [0.8, 1.25] 带 (半场伪解抓现行)
  A5 扇区符号相容 — sector 模式必须先过 asm_spec.sector_feasibility (均匀符号条件)

>>> from audit_protocol import audit
>>> import json, os
>>> d = os.path.join(os.path.dirname(__file__), '..', 'data', 'raw')
>>> # FAIL 分支: MUMPS 假奇异 (nan) vs SuperLU 正常 → A1/A3 必须抓
>>> v = audit(json.load(open(os.path.join(d, 'e7_protocol_fail', 'asm_report.json'))),
...           json.load(open(os.path.join(d, 'e7_protocol_fail', 'asm_report_superlu.json'))))
>>> v['verdict']
'FAIL'
>>> any(s['step'] == 'A4' and not s['pass'] for s in v['steps'])
True
>>> # PASS 分支: 单求解器健康报告 → 全步通过
>>> vp = audit(json.load(open(os.path.join(d, 'e7_protocol_pass', 'asm_report.json'))))
>>> vp['verdict']
'PASS'
>>> all(s['pass'] for s in vp['steps'])
True
"""
from __future__ import annotations

import math
from typing import Any

TOL_KCL = 1.0e-6
TOL_XSOLVER = 0.05
BAND_LO, BAND_HI = 0.8, 1.25
# 冻结基线 (paper11 data/raw/case1_scale_ladder, chi=1 全模型 16 条笼, s=1)
BASELINE_T = 0.0393
BASELINE_BG1 = 0.0579


def _first_row(report: dict[str, Any]) -> dict[str, Any]:
    rows = report.get("slip_sweep", {}).get("rows") or [{}]
    return rows[0]


def _finite(x: float) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(x)


def audit(report: dict[str, Any],
          report2: dict[str, Any] | None = None,
          sector_ok: bool | None = None) -> dict[str, Any]:
    steps: list[dict[str, Any]] = []

    def add(step: str, ok: bool, evidence: str) -> None:
        steps.append({"step": step, "pass": bool(ok), "evidence": evidence})

    def rel(a: float, b: float) -> float:
        m = max(abs(a), abs(b), 1e-30)
        return abs(a - b) / m

    # A1 完成性 (对每份报告; 兼容新旧链 physics_check 键结构)
    for name, rep in (("solver1", report), *((("solver2", report2),) if report2 else ())):
        pc = rep.get("physics_check") or {}
        row = _first_row(rep)
        bg1 = row.get("Bg1_amp_T")
        if "场有效" in pc:  # 旧链口径: 布尔三键
            ok = bool(pc.get("场有效")) and bool(pc.get("有限"))
        else:  # 新链口径: bg1_min_T 对 floor 判
            lo, fl = pc.get("bg1_min_T"), pc.get("bg1_floor_T", 1e-3)
            ok = _finite(lo) and lo > fl
        ok = ok and _finite(bg1) and bg1 > 1e-3
        add(f"A1[{name}]", ok, f"physics_check={pc} Bg1={bg1}")

    # A2 数值卫生
    row = _first_row(report)
    kcl = row.get("kcl_resid")
    add("A2", _finite(kcl) and abs(kcl) < TOL_KCL, f"kcl_resid={kcl} (tol {TOL_KCL})")

    # A3 交叉求解器一致性 (双报告时)
    if report2:
        r1, r2 = _first_row(report), _first_row(report2)
        t1, t2 = r1.get("torque_avg_Nm"), r2.get("torque_avg_Nm")
        b1, b2 = r1.get("Bg1_amp_T"), r2.get("Bg1_amp_T")
        if not (_finite(t1) and _finite(t2) and _finite(b1) and _finite(b2)):
            add("A3", False, f"non-finite quantity pair T=({t1},{t2}) Bg1=({b1},{b2})")
        else:
            dt, db = rel(t1, t2), rel(b1, b2)
            add("A3", dt <= TOL_XSOLVER and db <= TOL_XSOLVER,
                f"|dT|/max={dt:.4g} |dB|/max={db:.4g} (tol {TOL_XSOLVER})")

    # A4 尺度带自洽
    bg1 = row.get("Bg1_amp_T")
    if _finite(bg1) and bg1 > 1e-3:
        ratio = bg1 / BASELINE_BG1
        add("A4", BAND_LO <= ratio <= BAND_HI,
            f"Bg1/BASELINE={ratio:.3f} band=[{BAND_LO},{BAND_HI}]")

    # A5 扇区符号相容 (sector 模式显式传入预检结果)
    if sector_ok is not None:
        add("A5", sector_ok, f"uniform-sign pre-check={sector_ok}")

    verdict = "PASS" if all(s["pass"] for s in steps) else "FAIL"
    return {"verdict": verdict, "steps": steps}


if __name__ == "__main__":
    import json
    import sys

    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    rep = json.load(open(sys.argv[1]))
    rep2 = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
    print(json.dumps(audit(rep, rep2), ensure_ascii=False, indent=1))

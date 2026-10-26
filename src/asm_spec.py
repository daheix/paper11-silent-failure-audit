#!/usr/bin/env python3
"""ASM 笼型转差链口径卡 (doctest 可执行) — 数值口径唯一事实源

链路 (M-ASM-2 频域相量稳态): run_emag_asm.py →
      generate_geo.py(rotor_type='cage' → bar_k 导条区域) →
      motor_mag_asm.pro (GetDP Type Complex + Frequency = s×f_line, 固体导体
      全项耦合 + 端环电阻 Network) → asm_report.json

口径 (constants/asm_design.json 单一事实源):
- 转子系频域: 定子三相源相量频率 = s×f_line (转差频率), 转子网格静止 —
  每转差点单次线性解, 无时间推进 (天然消除 M-ASM-1 时域磁通俘获病)。
- 导条: 教程 8 固体导体全项耦合 (σDt{a}+σur/CoefGeo 双向+全局量), 频域
  DtDof 由 GetDP 自动替换 jω。
- 端环: 每端 n_bars 段集总电阻 R_seg = ρ·L_seg/A_ring, Network 环形电路
  (导条支路 {k, n+k}, 顶环 {k, k+1}, 底环 {n+k, n+k+1}, 节点 1..2n)。
- 频域线性铁磁: Flag_IronLinear=1 (Mur_Iron 定值, 未建谐波饱和, 披露)。
- 时均转矩 (相量=峰值约定): T = L·r²·∮ (1/2)Re(Bn·Bt*)/μ0 dθ。
- Bg1 幅值: 气隙圆周复数采样的 rfft 基波 |·|×2。
- 物理自检: T(s_max) ≥ ratio × T(s_min) (上升性) + Bg1 场量有效性下限。

>>> cfg = load_asm_spec()
>>> round(source_freq_hz(0.02, 50.0), 6)
1.0
>>> round(bar_area(3.0), 12) == round(3.141592653589793 * 1.5e-3**2, 12)
True
>>> round(bar_half_deg(3.0, 0.028), 6) > 0
True
>>> round(ring_resistance(3.5e7, 35.0, 4.0, 16), 7) > 0
True
>>> import math
>>> round(phasor_torque_density(1e-3, 1e-3) / (0.5e-6 / (4e-7 * math.pi)), 6)
1.0
>>> phasor_torque_density(complex(1, 0), complex(0, 2))
0.0
>>> round(time_avg_torque([1.0, -1.0, 1.0, -1.0], 0.1, 0.03), 9)
0.0
>>> round(time_avg_torque([2.0, 2.0], 0.1, 0.03), 9) == round(0.1 * 0.0009 * 2 * 6.283185307179586, 9)
True
>>> check_slip(0.5)
True
>>> check_slip(-0.2)
True
>>> check_slip(1.5)
Traceback (most recent call last):
    ...
ValueError: 转差 s=1.5 超出物理范围 [-0.5, 1]
>>> slip_mode(-0.2)
'generating (超同步发电, s<0)'
>>> slip_mode(0.05)
'motoring (亚同步电动, 0<s<1)'
>>> slip_mode(1.0)
'standstill 启动 (s=1)'
>>> progress_keys()
['stage', 'percent', 'elapsed_s', 'eta_s', 'msg', 'timestamp']
"""
import json
import math
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
# 打包自包含 (m20): doctest 在仓库源码树与捆绑包两处均可解析 constants;
# 环境变量 CHINASIM_MOTOR_CONSTANTS 可显式覆盖 (CI/部署场景)。
_CONST = os.environ.get("CHINASIM_MOTOR_CONSTANTS") or os.path.join(
    os.path.dirname(os.path.dirname(_HERE)), "constants")


def load_asm_spec():
    """ASM 口径默认值 — 单一事实源 constants/asm_design.json。"""
    path = os.path.join(_CONST, "asm_design.json")
    if not os.path.exists(path):
        # 捆绑包回退: 与本模块同目录的 asm_design.snapshot.json (发布物内打包)
        path = os.path.join(_HERE, "asm_design.snapshot.json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def source_freq_hz(slip, line_freq_hz):
    """转子系定子源相量频率 = slip × f_line (Hz)。"""
    return float(slip) * float(line_freq_hz)


def bar_cross_section_area_cage(bar_diameter_mm, r_ri_m, r_rm_m):
    """cage 路径导条实际截面积 [m²] = bar_half_deg·(R_rm²−R_ri²)/2 (扇形环带)。

    几何实现 (generate_geo.rotor_type='cage') 中 bar_k=转子带全厚扇形,
    bar_diameter_mm 为弧长当量直径 (只进角宽 arcsin(d/2/R_mid)) —
    M-ASM-10 复判: 3mm 薄带内 φ3 圆棒径向内切退化, 扇形=唯一合理实现,
    本函数为电路侧 (P_cuR 的 R') 与 FE 侧 (σ·A 体积) 的统一面积口径。

    >>> round(bar_cross_section_area_cage(3.0, 0.026, 0.029) * 1e6, 3)
    9.004
    """
    r_mid = 0.5 * (float(r_ri_m) + float(r_rm_m))
    half_deg = math.degrees(math.asin(
        (float(bar_diameter_mm) * 1e-3 / 2.0) / r_mid))
    return math.radians(2.0 * half_deg) * (float(r_rm_m) ** 2
                                           - float(r_ri_m) ** 2) / 2.0


def bar_area(bar_diameter_mm):
    """导条截面积 [m²] = π d²/4。"""
    return math.pi * (float(bar_diameter_mm) * 1e-3) ** 2 / 4.0


def bar_half_deg(bar_diameter_mm, r_mid_m):
    """导条角半宽 [°] = arcsin((d/2)/R_mid) — 几何生成边界角。"""
    return math.degrees(math.asin((float(bar_diameter_mm) * 1e-3 / 2.0)
                                  / float(r_mid_m)))


def ring_resistance(sigma_conductor, ring_radius_mm, ring_bar_diameter_mm,
                    n_bars):
    """端环单段电阻 [Ω] = ρ·L_seg/A_ring; L_seg = 2πR_ring/n_bars。"""
    rho = 1.0 / float(sigma_conductor)
    l_seg = 2.0 * math.pi * float(ring_radius_mm) * 1e-3 / float(n_bars)
    a_ring = bar_area(ring_bar_diameter_mm)
    return rho * l_seg / a_ring


def phasor_torque_density(bn, bt):
    """时均切向磁应力密度 [N/m²] = (1/2)Re(Bn·Bt*)/μ0 (相量=峰值约定)。

    bn/bt: 该点复数相量 (峰值)。峰值约定依据: 官方 machine_magstadyn_a.pro
    Torque_Maxwell_cplx 同含 0.5 且官方恒等式自洽 3.6% (R30 钟摆终审)。
    RMS 假说被否决: 磁路量级比 0.98 属磁路估算 ±30% 误差带内巧合
    (槽面积扇形近似/漏磁/磁压降均未计入)。"""
    mu0 = 4e-7 * math.pi
    return 0.5 * (complex(bn) * complex(bt).conjugate()).real / mu0


def time_avg_torque(p_tan_avg_list, l_stack_m, r_sample_m, scale=None):
    """Maxwell 应力时均转矩 [Nm] = L·r²·mean(p_tan)·2π。

    p_tan_avg_list: 采样圆逐点时均切向应力密度 (N/m²)。"""
    l = float(l_stack_m)
    r = float(r_sample_m)
    mean_p = sum(p_tan_avg_list) / len(p_tan_avg_list)
    return l * r * r * mean_p * 2.0 * math.pi


def ring_inductance(cfg):
    """端环每段集总电感 [H] — constants/asm_design.json l_endring_segment_H
    (官方 ASM 配方必需件: 纯阻环网数学退化; 官方锚 4.8nH)。"""
    return float(cfg["l_endring_segment_H"])


MU0 = 4e-7 * math.pi


def mu_eff_from_bmax(b_max_t, bh_table, mu_r_linear_cap=7000.0):
    """M-ASM-3 谐波饱和: 由局部磁密幅值查 B-H 表得有效相对磁导率。

    频域时均意义下的有效磁导率迭代法 (effective permeability iteration,
    行业标准工程口径): μ_eff = B_max / (μ0·H(B_max)), B-H 表线性内插。
    H 低于表首点 (线性区) 时用首点斜率并按 mu_r_linear_cap 截顶;
    B 超出表尾 (深度饱和) 时用尾点斜率外推。

    >>> mu_eff_from_bmax(0.05, [[0.0, 0.0], [0.9, 100.0]], 7000.0)
    7000.0
    >>> mu_eff_from_bmax(0.9, [[0.0, 0.0], [0.9, 100.0]], 7000.0)  # cap 生效
    7000.0
    >>> mu_eff_from_bmax(1.8, [[0.0, 0.0], [0.9, 100.0], [1.8, 4600.0]], 7000.0) > 0
    True
    """
    b = float(b_max_t)
    pts = sorted((float(p[0]), float(p[1])) for p in bh_table)
    if b <= pts[0][0]:
        h = pts[0][1] + (b - pts[0][0]) * _bh_slope(pts, 0)
    else:
        h = None
        for i in range(len(pts) - 1):
            b0, h0 = pts[i]
            b1, h1 = pts[i + 1]
            if b0 <= b <= b1:
                t = 0.0 if b1 == b0 else (b - b0) / (b1 - b0)
                h = h0 + t * (h1 - h0)
                break
        if h is None:
            h = pts[-1][1] + (b - pts[-1][0]) * _bh_slope(pts, len(pts) - 2)
    if h <= 0.0:
        return float(mu_r_linear_cap)
    return min(b / (MU0 * h), float(mu_r_linear_cap))


def _bh_slope(pts, i):
    """B-H 表第 i 段斜率 (H/B 视角: dH/dB)。"""
    db = pts[i + 1][0] - pts[i][0]
    dh = pts[i + 1][1] - pts[i][1]
    return dh / db if db > 0 else 0.0



def sector_feasibility(n_slots, n_poles, sector_deg, n_bars=None):
    """扇区化可行性判据 (paper11 均匀符号条件, M-ASM-17 机理终版)。

    周期约束 Link 的单一符号 (整周期 +1 / 反周期 -1) 必须对解中所有显著
    空间谐波阶 m 同时成立, 否则扇区化数学上不可行 (强制符号冲突 → 量值
    污染: 12s4p/90° 实测基波 x1.44/转矩 x1.54)。显著阶取基波 p、槽谐波
    k*n_slots±p 与齿谐波 k*n_slots (k=1..3), 笼条谐波 k*n_bars±p (若给 n_bars)。

    返回 dict(feasible, signs{m: 1|-1|None}, conflicts[...]);
    None = 非整非反 (直接不相容)。

    >>> f = sector_feasibility(12, 4, 90)
    >>> f["feasible"]
    False
    >>> f["signs"][2]     # 基波反周期
    -1
    >>> f["signs"][12]    # 齿谐波整周期, 与基波符号冲突
    1
    >>> sorted(f["conflicts"])[:3]
    [12, 24, 36]
    >>> sector_feasibility(12, 4, 180)["feasible"]   # 半圆: 符号一致
    True
    >>> sector_feasibility(12, 2, 180)["feasible"]   # 2 极: 基波-1 vs 齿谐+1 冲突
    False
    """
    p = n_poles // 2
    orders = {p}
    for k in (1, 2, 3):
        orders |= {k * n_slots - p, k * n_slots + p, k * n_slots}
    if n_bars:
        for k in (1, 2):
            orders |= {k * n_bars - p, k * n_bars + p}
    signs, conflicts = {}, []
    for m in sorted(orders):
        x = int(round(m * sector_deg)) % 360 if sector_deg else None
        s = {0: 1, 180: -1}.get(x)
        signs[m] = s
        if s is None or s != signs[p]:
            conflicts.append(m)
    return {"feasible": not conflicts, "signs": signs,
            "conflicts": sorted(set(conflicts))}


def check_slip(slip):
    """转差物理范围 [-0.5, 1] 校验 (s<0 = 超同步发电区, R19 扩展)。"""
    s = float(slip)
    if not (-0.5 <= s <= 1.0):
        raise ValueError(f"转差 s={s} 超出物理范围 [-0.5, 1]")
    return True


def slip_mode(slip):
    """转差象限标注 (电动/发电/启动) — 求解后报告归属用。

    >>> slip_mode(-0.2)
    'generating (超同步发电, s<0)'
    >>> slip_mode(0.05)
    'motoring (亚同步电动, 0<s<1)'
    >>> slip_mode(1.0)
    'standstill 启动 (s=1)'
    """
    s = float(slip)
    if s < 0.0:
        return "generating (超同步发电, s<0)"
    if s == 1.0:
        return "standstill 启动 (s=1)"
    return "motoring (亚同步电动, 0<s<1)"
    return True


def progress_keys():
    """_progress.json 键 schema (与静磁/瞬态链同, C++ 1s 轮询)。"""
    return ["stage", "percent", "elapsed_s", "eta_s", "msg", "timestamp"]


if __name__ == "__main__":
    import doctest
    failures, _ = doctest.testmod(verbose=False)
    raise SystemExit(1 if failures else 0)

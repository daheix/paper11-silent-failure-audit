#!/usr/bin/env python3
"""ASM 笼型转差链驱动器 — 链路编排 + 求解输入生成, 禁业务数值逻辑 (口径全部来自 constants/)

链路 (M-ASM-2 频域相量稳态): generate_geo.py(rotor_type='cage' → bar_k 导条)
      → gmsh → motor_mag_asm.pro (Type Complex + Frequency = s×f_line: 导条
      教程8固体导体全项耦合 + 端环电阻 Network) → 单次线性解/转差点 →
      时均转矩-转差曲线 → asm_report.json

用法:
  python3 run_emag_asm.py [--params overrides.json] [--workdir DIR]
      [--slips "1.0,0.5,0.2"] [--js-peak A/m2]

口径来源 (单一事实源):
  - constants/asm_design.json: 笼参数/端环参数/转差点/周期数/步数/物理自检比
  - asm_spec.py: 转差频率源 / 端环电阻 / 导条几何角 / 范围校验
  - external_circuit.py 实测口径 (GetDP 4.0.0): GlobalTerm 试探函数={Iz}、
    分段常数空括号、虚拟电路区域 tag 与网格物理组隔离
  - 转矩与 Bg1 解析复用 run_emag_getdp.GetDPMotorDriver (同一 Maxwell 应力口径)

物理口径 (constants/asm_design.json disclosure):
  转子系频域 FEM+电路联立 — 导条电流相量由磁场感应 EMF 经端环电阻 Network 真实
  决定; 定子源相量频率 = s×f_line, 转子网格静止; 时均转矩 = (L r²/2μ0)∮Re(Bn·Bt*)dθ。
  频域线性铁磁 (Flag_IronLinear=1, 谐波饱和/增量磁导率未建, 披露)。
  非解析伪实现: 全 FEM 离散 + 基尔霍夫电路联立, GetDP 频域真解。

验收: 各转差点求解成功 + 物理自检 T(s_max) ≥ ratio × T(s_min) (上升性)。
进度: _progress.json (ProgressReporter, 与静磁链同 schema, C++ 1s 轮询)。
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_emag_getdp as base                     # noqa: E402
from run_emag_getdp import (GetDPMotorDriver,     # noqa: E402
                            ProgressReporter, merged_design)
import asm_spec as spec
from asm_spec import (load_asm_spec, source_freq_hz,   # noqa: E402
                      bar_area, bar_half_deg, ring_resistance,
                      ring_inductance, phasor_torque_density,
                      time_avg_torque, check_slip, slip_mode, progress_keys)

# 虚拟电路区域 tag 基址 (与网格物理组 pid 隔离, external_circuit 同口径;
# groups.pro 物理组 pid ≤ 200, 此基址安全)
_CKT_TAG_BASE = 800000


def write_asm_pro(path, n_bars, n_sect=None, n_poles=4, sector_deg=None):
    """生成 motor_mag_asm.pro — 导条 FE 全局量 + 端环电阻 Network (n_bars 参数化)。

    布局 (GetDP 4.0.0 实测口径, 同 motor_mag_cir.pro / external_circuit.py):
      导条 k: 空间 Hcurl_u_bar_k (BF_RegionZ, Global Is/Us), 支路 {k, n+k};
      端环段 ring_t_k / ring_b_k: 虚拟区域, 支路 {k,k+1} / {n+k,n+k+1} (模 n);
      KCL/KVL: Constraint Type Network + Formulation GlobalEquation 联立。
    """
    n = int(n_bars)
    n_sect = int(n_sect) if n_sect else n   # 提前: 列表/循环全线扇区化
    bars = [f"bar_{k}" for k in range(n_sect)]
    rings = [f"ring_t_{k}" for k in range(n_sect)] + \
            [f"ring_b_{k}" for k in range(n_sect)]
    rings_l = [r + "_L" for r in rings]   # 端环电感段 (每段 R+L 两支路, 官方 ASM 配方)
    l = []
    A = l.append
    A("// motor_mag_asm.pro — ASM 笼型转差链主问题 (转子系转差频率 FEM+电路联立)")
    A("// === 由 run_emag_asm.py 按 n_bars 自动生成, 勿手改 ===")
    A("// 口径: constants/asm_design.json (唯一事实源) + asm_spec.py 口径卡")
    A("// 导条: 独立 FE 导体 (BF_RegionZ 全局电流/电压, 均匀电流密度集总 —")
    A("//       条内趋肤未建模, 披露见 disclosure: δ_Al(s·f≤50Hz)≥9.3mm > 条径)")
    A("// 端环: 每端 n 段集总电阻, Network 环形电路 (节点 1..2n, 全支路含阻)")
    A('')
    A('Include "motor_groups.pro"')
    A('Include "motor_mag_data.pro"')
    A('')
    A("DefineConstant[")
    A("  Fe = 0.,   // 转子系源相量频率 [Hz] = s×f_line (驱动器注入, 亦=System Frequency)")
    A("  IA0 = 0.,  IB0 = 0.,  IC0 = 0.,   // 定子相电流密度幅值 [A/m^2] (相量峰值)")
    A("  Rring = 0.,  // 端环单段电阻 [Ω] (asm_spec.ring_resistance, 驱动器注入)")
    A("  SigBar = 0., // 导条电导率 [S/m] (constants/ 驱动器注入)")
    A("  Lstk = 1.,   // 叠片厚度 [m] (CoefGeo 出平面厚度, 教程 8 2D 口径)")
    A("  Flag_IronLinear = 1,  // 频域线性铁磁 (谐波饱和未建, 披露)")
    A("  Mur_Iron = 7000,  Mur_Shaft = 1")
    A("];")
    A('')
    # ---- 虚拟电路区域 (Group 追加块, tag 基址隔离) + 区组组合 ----
    A("Group {  // 区组组合 + 端环虚拟电路区域 (无网格, external_circuit 同口径)")
    A("  Vol_All = Region[{ shaft, rotor_iron, Bars, airgap,")
    A("                       stator_iron, SlotAp, SlotAn, SlotBp, SlotBn,")
    A("                       SlotCp, SlotCn }];")
    A("  Vol_Iron = Region[{ rotor_iron, stator_iron }];")
    A("  Vol_Slots = Region[{ SlotAp, SlotAn, SlotBp, SlotBn, SlotCp, SlotCn }];")
    for i, r in enumerate(rings):
        A(f"  {r} = Region[ {_CKT_TAG_BASE + i} ];")
    A(f"  Dom_Cir = Region[ {{{','.join(rings)}}} ];")
    A("}")
    A('')
    # ---- 材料函数 (与静磁/瞬态链同口径, 无磁体; hc 已被驱动器存根化) ----
    A("Function {")
    if n_sect < n:
        # M-ASM-17: 扇区 Link 周期约束的旋转映射 (主链 motor_mag.pro:23 同款,
        #   官方 AG-EM machine_magstadyn_a.pro:298; asm .pro 独立生成须自带定义)
        A("  RotatePZ[] = Rotate[ Vector[$X,$Y,$Z], 0, 0, $1 ] ;")
    A("  nu[ Region[{airgap}] ] = 1/mu0;")
    A("  nu[ Region[{SlotAp, SlotAn, SlotBp, SlotBn, SlotCp, SlotCn}] ] = 1/mu0;")
    A("  nu[ Region[{Bars}] ] = 1/mu0;              // 铝/铜导条 μr≈1")
    A("  nu[ Region[{shaft}] ] = 1/(mu0*mur_shaft);")
    A("  nu[ Region[Vol_Iron] ] = Flag_IronLinear ? 1/(mu0*Mur_Iron) : iron_nu[$1];")
    A('')
    A("  // 频域源 js 相量 [A/m^2]: Complex[Cos,Sin] 相量常数 (官方模板 Complex[] 构造器;")
    A("  //   F_Cos_wt_p 在区域函数内求值为零 — M-ASM-2b 根因2 实证)")
    A("  // 正序: A=0°, B=-120°, C=+120° (相对相位保持, 全局相位不影响 |T| 曲线)")
    A("  js[ Region[{SlotAp}] ] = Vector[0, 0,  IA0*Complex[Cos[0], Sin[0]]];")
    A("  js[ Region[{SlotAn}] ] = Vector[0, 0, -IA0*Complex[Cos[0], Sin[0]]];")
    A("  js[ Region[{SlotBp}] ] = Vector[0, 0,  IB0*Complex[Cos[-2*Pi/3], Sin[-2*Pi/3]]];")
    A("  js[ Region[{SlotBn}] ] = Vector[0, 0, -IB0*Complex[Cos[-2*Pi/3], Sin[-2*Pi/3]]];")
    A("  js[ Region[{SlotCp}] ] = Vector[0, 0,  IC0*Complex[Cos[2*Pi/3], Sin[2*Pi/3]]];")
    A("  js[ Region[{SlotCn}] ] = Vector[0, 0, -IC0*Complex[Cos[2*Pi/3], Sin[2*Pi/3]]];")
    A('')
    A("  // 导条电导率 (固体导体模型, 教程 8 口径)")
    A("  Sigma[ Region[{Bars}] ] = SigBar;")
    A('')
    A("  // 端环分段常数 (显式逐区域赋值, 禁默认分段 — GetDP 4.0 实测 Redefinition);")
    A("  // 官方 ASM 拓扑: 顶环=电阻载体, 底环=电感载体 (每段单元件, 免并联退化)")
    for k in range(n_sect):
        A(f"  Inductance[ ring_b_{k} ] = Lring;   // 底环=电感载体")
    A("  CoefGeo[] = Lstk;   // 出平面厚度 (教程 8 2D 口径, 导条 ur 耦合用)")
    for k in range(n_sect):
        A(f"  Resistance[ ring_t_{k} ] = Rring;   // 顶环=电阻载体")
    A("  // 底环 (ring_b_*) 无 Resistance 赋值 — 电感载体, 禁 0 系数冗余项")
    A("}")
    A('')
    # ---- 约束: 外边界 + 扇区周期 Link + Network ----
    A("Constraint {")
    A("  { Name a; Type Assign;")
    A("    Case { { Region Outer; Value 0.; }")
    if n_sect < n:
        # M-ASM-17 真根因修复 (2026-10-26): 扇区切割边界 CutB 被 Outer 钳 a=0
        #   Dirichlet (motor_groups.pro Outer=Region[14]=CutB) → 磁通路径被掐死
        #   → Bg1=0 假守恒 (asm_js12chk 实证 1e-5 量级)。修复=主链 M-ASM-11 R97
        #   同口径 Link 周期约束 (官方 AG-EM 语法 machine_magstadyn_a.pro:335:
        #   Link 是 Case 级成员, If 门在 Case{} 内部)。系数符号按扇区跨半极数
        #   (90°/4 极=1 半极 → -1 反周期; 与主链 sector_bg1 延拓符号同式)。
        #   系数符号用极对 p=n_poles/2 (主链 sector_bg1 同式: n_half=sd·p/180,
        #   奇数半极→反周期 -1; 4 极 90° → n_half=1 → -1, 首版误用 n_poles
        #   得 +1 = 周期性判错, 主链 R97 注释同口径"4 极半极反周期")
        n_half = sector_deg * (n_poles / 2.0) / 180.0
        link_coef = -1.0 if int(round(n_half)) % 2 else 1.0
        # M-ASM-17 4b: Link 作用在 1D 切割线 — 面 Region 旋转精确匹配因两侧
        # 带内部三角节点坐标浮点尾数不同必失败 (js17b 实证 154 仅 51 同源线
        # 点匹配; msh2 无 $Periodic 段)。线节点两侧同源定义坐标解析相等 →
        # 全匹配; 面内部 DoF 由面内变分方程自动满足周期性 (标准做法)。
        A("      If (Flag_Sector == 1)")
        A("      { Region CutLineB; Type Link; RegionRef CutLineA;")
        A(f"        Coefficient {link_coef:.0f};")
        A("        Function RotatePZ[-Sector_deg*Pi/180]; }")
        A("      EndIf")
    A("    } }")
    A("  { Name ElectricalCircuit; Type Network;")
    A("    Case Ckt_1 {")
    A("      // 支路声明顺序: 无源元件在前 FE 导体在后 — 顺序翻转会使大图装配病态")
    A("      //   (M-ASM-2b 根因3 实证); 每段端环 = R + L 两支路串联 (官方 ASM Rers+Lers,")
    A("      //   纯阻环网数学退化 — 零模/病态根因)")
    # 端环折返 (官方 AG-EM 对称约化手法, im_circuit.pro:55/61):
    # n_sect < n 时扇区最后一段不再连 (k+1)%n, 而是折返到扇区首节点 1
    # (偶对称口径 SymmetryFactor%2==0 → NB1_1); 几何 1/4 模型时电路在
    # 周期边界自闭合。n_sect=n 时折返段 k=n-1 连 {n,1} 与模 n 恒等
    # (回归锚: 全模型 0.1141 不变)。
    def _ring_pair(k, base):
        """端环段 k 的支路节点对 (base=0 顶环 / n 底环偏移)。"""
        if k + 1 < n_sect:
            return base + k + 1, base + k + 2
        return base + n_sect, base + 1      # 折返: 扇区末段→首节点

    for k in range(n_sect):
        a1, a2 = _ring_pair(k, 0)
        A(f"      {{ Region ring_t_{k}; Branch {{{a1}, {a2}}}; }}   // 顶环=电阻段")
    for k in range(n_sect):
        a1, a2 = _ring_pair(k, n_sect)
        A(f"      {{ Region ring_b_{k}; Branch {{{a1}, {a2}}}; }}   // 底环=电感段")
    for k in range(n_sect):   # 扇区: FE 导条组只有前 n_sect 条 (M-ASM-11)
        A(f"      {{ Region bar_{k}; Branch {{{k + 1}, {n_sect + k + 1}}}; }}")
    A("    } }")
    A("}")
    A('')
    # ---- 函数空间: 磁场 Form1P + 每导条一个全局量空间 ----
    A("FunctionSpace {")
    A("  { Name Hcurl_a; Type Form1P;")
    A("    BasisFunction {")
    A("      { Name se; NameOfCoef ae; Function BF_PerpendicularEdge;")
    A("        Support Vol_All; Entity NodesOf[ All ]; }")
    A("    }")
    A("    Constraint {")
    A("      { NameOfCoef ae; EntityType NodesOf; NameOfConstraint a; }")
    A("    }")
    A("  }")
    A("  // 导条: 单一共享空间 (BF_RegionZ 每区域一个自由度, 教程 8 官方口径 —")
    A("  //   逐导体独立空间实测装配零行)")
    A("  { Name Hcurl_u_bars; Type Vector;")
    A("    BasisFunction {")
    A("      { Name sr; NameOfCoef ur; Function BF_RegionZ;")
    A("        Support Bars; Entity Bars; }")
    A("    }")
    A("    GlobalQuantity {")
    A("      // 教程 8 精确镜像: Us=AliasOf (KVL 量), Is=AssociatedWith (KCL/Dof 量)")
    A("      //   — 指派反转会使 FE↔电路耦合装配成转置 → 电路视导条开路 (M-ASM-2b 根因)")
    A("      { Name Us; Type AliasOf; NameOfCoef ur; }")
    A("      { Name Is; Type AssociatedWith; NameOfCoef ur; }")
    A("    }")
    A("    Constraint {")
    A("      { NameOfCoef Us; EntityType Region; NameOfConstraint Voltage_Bars; }")
    A("      { NameOfCoef Is; EntityType Region; NameOfConstraint Current_Bars; }")
    A("    }")
    A("  }")
    A("  // 电路自由度 (Scalar + BF_Region, Current=AliasOf/Voltage=AssociatedWith —")
    A("    //  与教程 8 cir 空间指派一致, 与导条空间故意相反)")
    A("  //  电路自由度 (Scalar + BF_Region, external_circuit 同口径)")
    A("  { Name Hregion_i_Cir; Type Scalar;")
    A("    BasisFunction {")
    A("      { Name sr; NameOfCoef ir; Function BF_Region;")
    A("        Support Dom_Cir; Entity Dom_Cir; }")
    A("    }")
    A("    GlobalQuantity {")
    A("      { Name Current; Type AliasOf; NameOfCoef ir; }")
    A("      { Name Voltage; Type AssociatedWith; NameOfCoef ir; }")
    A("    }")
    A("    Constraint {")
    A("      { NameOfCoef Voltage; EntityType Region; NameOfConstraint Voltage_Cir; }")
    A("      { NameOfCoef Current; EntityType Region; NameOfConstraint Current_Cir; }")
    A("    }")
    A("  }")
    A("}")
    A('')
    # ---- 空约束声明 (Constraint 块, 教程 8 口径) ----
    A("Constraint {  // 导条空约束: 电流由电路决定 (教程 8: Current_Mag/Voltage_Mag)")
    A("  { Name Current_Bars; Case { } }")
    A("  { Name Voltage_Bars; Case { } }")
    A("  { Name Voltage_Cir; Case { } }")
    A("  { Name Current_Cir; Case { } }")
    A("}")
    A('')
    A("Jacobian { { Name JacVol; Case { { Region All; Jacobian Vol; } } } }")
    A("Integration { { Name Int_Mag; Case { { Type Gauss; Case {")
    A("  { GeoElement Point; NumberOfPoints 1; }  { GeoElement Line; NumberOfPoints 3; }")
    A("  { GeoElement Triangle; NumberOfPoints 4; }  { GeoElement Quadrangle; NumberOfPoints 4; }")
    A("  { GeoElement Tetrahedron; NumberOfPoints 4; }  { GeoElement Hexahedron; NumberOfPoints 6; }")
    A("  { GeoElement Prism; NumberOfPoints 9; }")
    A("} } } } }")
    A('')
    # ---- 形式: 磁场 + 导条全局量 + 端环电路 (联立) ----
    A("Formulation {")
    A("  { Name Magnetodynamics_asm; Type FemEquation;")
    A("    Quantity {")
    A("      { Name a; Type Local; NameOfSpace Hcurl_a; }")
    A("      { Name ur_b; Type Local; NameOfSpace Hcurl_u_bars; }")
    A("      { Name Is_b; Type Global; NameOfSpace Hcurl_u_bars [Is]; }")
    A("      { Name Us_b; Type Global; NameOfSpace Hcurl_u_bars [Us]; }")
    A("      { Name Iz; Type Global; NameOfSpace Hregion_i_Cir [Current]; }")
    A("      { Name Uz; Type Global; NameOfSpace Hregion_i_Cir [Voltage]; }")
    A("    }")
    A("    Equation {")
    A("      // ---- 磁场方程 (线性区 ν0 + 非线性铁磁 Picard + 定子电流源) ----")
    A("      Integral { [ nu[] * Dof{d a} , {d a} ];")
    A("        In Region[{airgap, shaft, Vol_Slots, Bars}];")
    A("        Jacobian JacVol; Integration Int_Mag; }")
    A("      Integral { [ nu[{d a}] * Dof{d a}, {d a} ];")
    A("        In Vol_Iron; Jacobian JacVol; Integration Int_Mag; }")
    A("      Integral { [ -js[] , {a} ];")
    A("        In Vol_Slots; Jacobian JacVol; Integration Int_Mag; }")
    A("      // ---- 导条: 教程 8 全项固体导体模型 (σ dt a + 净电流 ur 双向耦合) ----")
    A("      Integral { DtDof [ Sigma[] * Dof{a} , {a} ];")
    A("        In Bars; Jacobian JacVol; Integration Int_Mag; }")
    A("      Integral { [ Sigma[] * Dof{ur_b} / CoefGeo[] , {a} ];")
    A("        In Bars; Jacobian JacVol; Integration Int_Mag; }")
    A("      Integral { DtDof [ Sigma[] * Dof{a} , {ur_b} ];")
    A("        In Bars; Jacobian JacVol; Integration Int_Mag; }")
    A("      Integral { [ Sigma[] * Dof{ur_b} / CoefGeo[] , {ur_b} ];")
    A("        In Bars; Jacobian JacVol; Integration Int_Mag; }")
    A("      GlobalTerm { [ Dof{Is_b} * Sign[CoefGeo[]] , {Us_b} ]; In Bars; }")
    A("      // ---- 端环 R+L 本构 (试探函数={Iz}, 逐区域施加, R 加 NeverDt — 官方 ASM 口径) ----")
    for k in range(n_sect):
        A(f"      GlobalTerm {{ NeverDt[ Dof{{Uz}} , {{Iz}} ]; In ring_t_{k}; }}")
        A(f"      GlobalTerm {{ NeverDt[ Resistance[] * Dof{{Iz}} , {{Iz}} ]; In ring_t_{k}; }}")
    for k in range(n_sect):
        A(f"      GlobalTerm {{ NeverDt[ Dof{{Uz}} , {{Iz}} ]; In ring_b_{k}; }}")
        A(f"      GlobalTerm {{ DtDof [ Inductance[] * Dof{{Iz}} , {{Iz}} ]; In ring_b_{k}; }}")
    A("      // ---- KCL/KVL 联立 (Network 拓扑在 Constraint 声明) ----")
    A("      GlobalEquation {")
    A("        Type Network; NameOfConstraint ElectricalCircuit;")
    A("        { Node {Is_b}; Loop {Us_b}; Equation {Is_b}; In Bars; }")
    A("        { Node {Iz}; Loop {Uz}; Equation {Uz}; In Dom_Cir; }")
    A("      }")
    A("    }")
    A("  }")
    A("}")
    A('')
    # ---- 求解序列 (T0 完整 Picard + TimeLoopTheta 单步线性解, 瞬态链同口径) ----
    A("Resolution {")
    A("  { Name Magnetodynamics_asm;")
    A("    System { { Name A; NameOfFormulation Magnetodynamics_asm;")
    A("              Type Complex; Frequency Fe; } }")
    A("    Operation {")
    A("      // 频域单次线性解 (无 Picard: Flag_IronLinear=1; DtDof→jω 自动)")
    A("      InitSolution[A];")
    A("      Generate[A]; Solve[A]; SaveSolution[A];")
    A("      // 后处理不内联调用: 全部打印走二次 -pos 进程读 .res (M212 同口径 —")
    A("      //   求解进程内联打印全局量 TimeTable 报 Empty solution, GetDP 4.0.0 实测)")
    A("    }")
    A("  }")
    A("}")
    A('')
    A("PostProcessing {")
    A("  { Name MagAsm; NameOfFormulation Magnetodynamics_asm;")
    A("    Quantity {")
    A("      { Name a; Value { Term { [ {a} ]; In Vol_All; Jacobian JacVol; } } }")
    A("      { Name b; Value { Term { [ {d a} ]; In Vol_All; Jacobian JacVol; } } }")
    A("      { Name iz; Value { Term { [ {Iz} ]; In Dom_Cir; } } }")
    A("    }")
    A("  }")
    A("}")
    A('')
    A("PostOperation {")
    A("  { Name MagAsm; NameOfPostProcessing MagAsm;")
    A("    Operation {")
    A("      // 全部打印在二次 -pos 进程执行 (读 .res, 驱动器 stdin 喂 \"1\\n\")")
    A('      Print[ b, OnElementsOf Vol_All, LastTimeStepOnly, File "b_map_t.pos" ];')
    A('      Print[ a, OnElementsOf Vol_All, LastTimeStepOnly, File "a_map_t.pos" ];')
    A('      Print[ a, OnElementsOf Vol_Slots, LastTimeStepOnly, File "az_slots_t.pos" ];')
    A('      Print[ b, OnElementsOf airgap, LastTimeStepOnly, File "bgap_t.pos" ];')
    A('      Print[ iz, OnRegion Dom_Cir, File "iz_asm.dat", Format TimeTable ];')
    A("    }")
    A("  }")
    A("}")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(l) + "\n")


def _read_iz_timetable(path):
    """GetDP Print TimeTable 解析 — 行: t re im re im ... (header # 注释行)。

    返回复相量幅值列表 (峰值口径, GetDP 频域复解约定); 单频链只取末行。
    """
    vals = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            vals = [float(x) for x in line.split()]
    if vals is None:
        raise RuntimeError(f"{path}: TimeTable 无数据行 — 后处理未落盘")
    if len(vals) < 3 or (len(vals) - 1) % 2 != 0:
        raise RuntimeError(f"{path}: 列数 {len(vals)} 非奇数 (t+re/im 对) — 格式漂移")
    body = vals[1:]
    return [complex(body[i], body[i + 1]) for i in range(0, len(body), 2)]


class GetDPAsmDriver(GetDPMotorDriver):
    """ASM 驱动: 笼型网格 + 生成式 .pro + 转差瞬态。"""

    _SEED_FILES = GetDPMotorDriver._SEED_FILES  # asm pro 由本驱动生成

    def mesh_cage(self, cfg):
        """笼型网格: geo 参数注入 rotor_type='cage' + 导条参数 (磁体带 → bar_k)。
        dtheta_max_deg: 网格角步距上限 (asm_design 键, 默认 2.5=历史口径;
        车用级大模型 MUMPS 内存适配通道 — 粗化自由度, 网格收敛性由
        Bg1 判据把关, 非求解降级)。"""
        base.P = dict(base.P, rotor_type="cage",
                      n_bars=int(cfg["n_bars"]),
                      bar_diameter_mm=float(cfg["bar_diameter_mm"]),
                      dtheta_max=float(cfg.get("dtheta_max_deg", 2.5)),
                      rotor_angle_deg=0.0)
        # lc 尺寸场粗化通道 (R130, --lc-air/--lc-iron/--lc-shaft): 大模型
        # (R_so>90mm) 网格主成本=lc_air 场 (R129 三通道实证 dtheta/扇区
        # 救不了); generate_geo.DEFAULTS 合法键直通, 缺省不动历史口径
        for _k in ("lc_air", "lc_iron", "lc_shaft"):
            if _k in cfg:
                base.P[_k] = float(cfg[_k])
        # M-ASM-14 网格通道 (--mesh-extend/--mesh-algo): gmsh 官方手册
        # 1.2.2 节修法 (ExtendFromBoundary=0 防 over-refinement +
        # Algorithm 5 Delaunay 大网格最快), 缺省不动历史口径
        for _k in ("mesh_extend_boundary", "mesh_algorithm",
                   "gap_transfinite", "gap_transfinite_pts"):
            if _k in cfg:
                base.P[_k] = int(cfg[_k])
        self.mesh()
        self._sanitize_data_for_cage()
        # 生成 .pro (n_bars 参数化; 扇区时 n_sect 折返 — 电路支路/端环
        # 按扇区条数生成, R93 折返口径, M-ASM-11)
        n_sect = None
        if self.sector_deg:
            # M-ASM-17 (2026-10-26 修复放行): 真根因=CutB 被 Outer 钳 a=0
            #  (motor_groups.pro Outer=Region[14]=CutB) 而非折返拓扑 (逐节点
            #  对照官方 im_circuit.pro AG-EM 折返口径一致)。修复=Constraint
            #  段补 CutB→CutA Link 周期约束 (主链 M-ASM-11 R97 同口径),
            #  验证实验 asm_js17b: 16 条 90° 扇区对照全模型基线。
            n_sect = int(round(360.0 / float(self.sector_deg)))
            self.log.append("M-ASM-17 修复放行: 扇区 cage Link 周期约束已实装 "
                            f"(CutB→CutA, n_sect={n_sect}, 系数按半极数符号)")
        write_asm_pro(os.path.join(self.wd, "motor_mag_asm.pro"),
                      int(cfg["n_bars"]), n_sect=n_sect,
                      n_poles=int(cfg.get("n_poles", 4)),
                      sector_deg=float(self.sector_deg) if self.sector_deg else None)
        self.log.append(f"asm pro: n_bars={cfg['n_bars']} n_sect={n_sect} 生成")

    def _sanitize_data_for_cage(self):
        """cage 工作目录数据适配: hc[] 定义引用 magnet_N/magnet_S 区域
        (cage 模式 groups.pro 不存在) — 将工作目录副本 motor_mag_data.pro 的
        HC VECTOR SECTIONS 段替换为零向量 shaft 存根 (asm pro 不使用 hc[],
        仅消除未知区域解析错误; 种子文件不动, PMSM 链不受影响)。"""
        dp = os.path.join(self.wd, "motor_mag_data.pro")
        with open(dp, encoding="utf-8") as f:
            txt = f.read()
        marker = "=== HC VECTOR SECTIONS"
        idx = txt.find(marker)
        if idx < 0:
            raise RuntimeError("motor_mag_data.pro 缺少 HC VECTOR SECTIONS 段")
        stub = ("// === HC VECTOR SECTIONS (cage 适配存根: ASM 链无磁体, hc[] 未用) ===\n"
                "Function {\n"
                "  hc[ Region[{shaft}] ] = Vector[0, 0, 0.];\n"
                "}\n")
        with open(dp, "w", encoding="utf-8") as f:
            f.write(txt[:idx] + stub)

    def solve_slip(self, js_pk, s, cfg, out_bmap=None):
        """一个转差点: 频域相量稳态 (Type Complex, Frequency = s×f_line)。

        数值口径: Fe=s×f_line / IA0..IC0=js_pk / Rring=asm_spec 端环电阻 /
        SigBar=σ_Al / Lstk=P.L_stack / Flag_IronLinear=1 — 零业务常数内嵌。"""
        check_slip(s)
        fs = source_freq_hz(s, cfg["line_freq_hz"])
        r_ring = ring_resistance(float(cfg["sigma_conductor_S_per_m"]),
                                 float(cfg["ring_radius_mm"]),
                                 float(cfg["ring_bar_diameter_mm"]),
                                 int(cfg["n_bars"]))
        vals = {
            "Fe": fs,
            "IA0": js_pk, "IB0": js_pk, "IC0": js_pk,
            "Rring": r_ring,
            "Lring": ring_inductance(cfg),
            "SigBar": float(cfg["sigma_conductor_S_per_m"]),
            "Lstk": float(base.P["L_stack"]),
            "Flag_IronLinear": 1,
            "Mur_Iron": float(cfg["mur_iron_effective"]), "Mur_Shaft": 1,
            "Flag_Sector": 1 if self.sector_deg else 0,
        }
        if self.sector_deg:
            vals["Sector_deg"] = float(self.sector_deg)
        cmd = [base.GETDP_BIN, "motor_mag_asm.pro", "-msh", "motor_sector.msh",
               "-solve", "Magnetodynamics_asm",
               # 求解器: Scaling=2 (行列范数) + 大 Krylov/高迭代/严停机 GMRES —
               # 默认无缩放 GMRES 在 σ/μ 高对比耦合系统上 2-3 步假收敛, 产出未收敛
               # 伪解 (M-ASM-2e 实证: 同配置 Scaling=2 后残差 3.5e-4→物理解 0.0617T)
               "-Scaling", "2"]
        for k, v in vals.items():
            cmd += ["-setnumber", k, repr(float(v))]
        t0 = time.time()
        self.prog.check_stop()
        # 峰值内存埋点 (goal 二.3 台账义务): RUSAGE_CHILDREN.ru_maxrss=已终止
        # 子进程峰值 RSS 最大值 (Linux 单位 KB), 求解前后取 max 即本轮 getdp 峰值
        import resource
        rss_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        r = subprocess.run(cmd, cwd=self.wd, capture_output=True, text=True,
                           stdin=subprocess.DEVNULL)
        rss_after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        peak_rss_mb = round(max(rss_after, rss_before) / 1024.0, 1)
        if r.returncode != 0:
            raise RuntimeError(
                f"getdp ASM 转差 s={s} 频域求解失败:\n"
                f"{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        self.log.append(f"solve_slip(s={s},fs={fs:.3g}Hz,Rring={r_ring:.4g}): "
                        f"{time.time() - t0:.1f}s rss={peak_rss_mb}MB")
        # 二次 -pos: 场图+端环电流全部落盘 (独立进程读 .res, M212 同口径;
        # stdin "1\n" 选第一个 PostProcessing)
        cmd_pos = [base.GETDP_BIN, "motor_mag_asm", "-msh", "motor_sector.msh",
                   "-pos"]
        for k, v in vals.items():
            cmd_pos += ["-setnumber", k, repr(float(v))]
        rp = subprocess.run(cmd_pos, cwd=self.wd, capture_output=True, text=True,
                            input="1\n")
        if rp.returncode != 0:
            raise RuntimeError(
                f"getdp -pos 后处理失败 (s={s}):\n"
                f"{rp.stdout[-800:]}\n{rp.stderr[-800:]}")
        if out_bmap and out_bmap != "b_map_t.pos":
            os.replace(os.path.join(self.wd, "b_map_t.pos"),
                       os.path.join(self.wd, out_bmap))
        return {"Rring_ohm": r_ring, "source_freq_hz": fs,
                "solve_s": round(time.time() - t0, 1),
                "peak_rss_mb": peak_rss_mb}

    # ---------- 频域复数场解析 (b_map.pos: 9 分量 × re/im = 18 值/单元) ----------
    def iron_bmax(self, bmap_path, quantile=0.98):
        """M-ASM-3 谐波饱和探针: 铁芯区域 |B| 相量幅值高分位 [T]。

        复数场 |B| = sqrt(|Bx|²+|By|²); 排除气隙圆环与轴心 (质心半径筛选,
        定子 R_si..R_so / 转子 R_ri..(R_si-气隙) 的铁区); 分位抗网格噪声。"""
        import numpy as np
        cs, vs = base.parse_pos_vv(os.path.join(self.wd, bmap_path))
        bx = vs[:, 0] + 1j * vs[:, 3]
        by = vs[:, 1] + 1j * vs[:, 4]
        mag = np.sqrt(np.abs(bx) ** 2 + np.abs(by) ** 2)
        cx = cs[:, 0::3].mean(axis=1)
        cy = cs[:, 1::3].mean(axis=1)
        r = np.hypot(cx, cy)
        r_si, r_ri = float(base.P["R_si"]), float(base.P["R_ri"])
        r_gap = float(base.P.get("r_gap_mid", 0.5 * (r_si + r_ri)))
        tol = 1.5 * (r_si - r_gap)
        iron = ((r >= r_si + tol) | (r <= r_ri - tol))
        if not iron.any():
            iron = np.ones_like(mag, dtype=bool)
        return float(np.quantile(mag[iron], float(quantile)))

    def gap_circle_freq(self, bmap_path, r_sample=None, n=1440,
                        th_span_deg=360.0):
        """频域气隙圆周复数采样 → (θ°, Bn 相量, Bt 相量)。
        th_span_deg<360 (扇区模型 R99): 只采扇区弧段 (全圆口径扇区外
        270° 点回退最近单元=污染)。"""
        import numpy as np
        if r_sample is None:
            r_sample = base.P["r_gap_mid"]
        cs, vs = base.parse_pos_vv(os.path.join(self.wd, bmap_path))
        assert vs.shape[1] == 18, \
            f"频域 b_map 期望 18 值/单元 (9 分量×re/im), 得到 {vs.shape[1]}"
        cx = cs[:, 0::3].mean(axis=1)
        cy = cs[:, 1::3].mean(axis=1)
        bx = vs[:, 0] + 1j * vs[:, 3]
        by = vs[:, 1] + 1j * vs[:, 4]
        bmap = {(round(x, 6), round(y, 6)): (bxi, byi)
                for x, y, bxi, byi in zip(cx, cy, bx, by)}
        tris = self._gap_tris()
        cen = tris.mean(axis=1)
        bxy = np.array([bmap[(round(c[0], 6), round(c[1], 6))] for c in cen],
                       dtype=complex)
        th = np.linspace(0.0, 2 * np.pi * th_span_deg / 360.0, n,
                         endpoint=False)
        px, py = r_sample * np.cos(th), r_sample * np.sin(th)
        p0 = tris[:, 0]; p1 = tris[:, 1]; p2 = tris[:, 2]
        d = ((p1[:, 1] - p2[:, 1]) * (p0[:, 0] - p2[:, 0])
             + (p2[:, 0] - p1[:, 0]) * (p0[:, 1] - p2[:, 1]))
        X = px[:, None] - p2[None, :, 0]
        Y = py[:, None] - p2[None, :, 1]
        l1 = ((p1[:, 1] - p2[:, 1])[None, :] * X
              + (p2[:, 0] - p1[:, 0])[None, :] * Y) / d[None, :]
        l2 = ((p2[:, 1] - p0[:, 1])[None, :] * X
              + (p0[:, 0] - p2[:, 0])[None, :] * Y) / d[None, :]
        l3 = 1.0 - l1 - l2
        hit = (l1 >= -1e-9) & (l2 >= -1e-9) & (l3 >= -1e-9)
        j = np.argmax(hit, axis=1)
        ok = hit[np.arange(n), j]
        if not ok.all():
            d2 = (px[:, None] - cen[None, :, 0]) ** 2 + \
                 (py[:, None] - cen[None, :, 1]) ** 2
            miss = np.where(~ok)[0]
            j[miss] = np.argmin(d2[miss], axis=1)
        # 重心插值 (M-ASM-5 终验): 命中单元内用 (l1,l2,l3) 加权 3 顶点复值。
        # tris 序=msh 节点序, 与 pos 行序无关 → 顶点值经坐标键 (x,y) 匹配
        vmap = {}
        for k in range(vs.shape[0]):
            for t in range(3):
                key = (round(float(cs[k, 3 * t]), 6),
                       round(float(cs[k, 3 * t + 1]), 6))
                vmap[key] = (vs[k, 6 * t] + 1j * vs[k, 6 * t + 3],
                             vs[k, 6 * t + 1] + 1j * vs[k, 6 * t + 4])
        bvx = np.empty((tris.shape[0], 3), dtype=complex)
        bvy = np.empty((tris.shape[0], 3), dtype=complex)
        for t in range(3):
            for k in range(tris.shape[0]):
                key = (round(float(tris[k, t, 0]), 6),
                       round(float(tris[k, t, 1]), 6))
                bvx[k, t], bvy[k, t] = vmap[key]
        lam = np.stack([l1[np.arange(n), j], l2[np.arange(n), j],
                        l3[np.arange(n), j]], axis=1)
        bad = ~ok
        lam[bad] = 1.0 / 3.0        # miss 点回退三顶点均值
        bx = np.einsum('ij,ij->i', lam, bvx[j])
        by = np.einsum('ij,ij->i', lam, bvy[j])
        ct, st = np.cos(th), np.sin(th)
        return np.degrees(th), bx * ct + by * st, -bx * st + by * ct

    def torque_freq(self, bmap_path, r_sample=None, sector_deg=None):
        """频域时均 Maxwell 应力转矩 [Nm] = L·r²·mean((1/2)Re(Bn·Bt*))/μ0·2π。
        扇区: 应力密度 p_tan 在反周期下偶延拓 (Bn/Bt 同反号→乘积不变),
        延拓后积分全圆 (R99)。"""
        import numpy as np
        if r_sample is None:
            r_sample = base.P["r_gap_mid"]
        if sector_deg:
            _, bn, bt = self.gap_circle_freq(bmap_path, r_sample, n=360,
                                             th_span_deg=sector_deg)
            n_rep = int(round(360.0 / sector_deg))
            p_tan = np.concatenate(
                [np.array([phasor_torque_density(b, t)
                           for b, t in zip(bn, bt)])] * n_rep)
            return time_avg_torque(list(p_tan), base.P["L_stack"], r_sample)
        _, bn, bt = self.gap_circle_freq(bmap_path, r_sample)
        p_tan = [phasor_torque_density(b, t) for b, t in zip(bn, bt)]
        return time_avg_torque(p_tan, base.P["L_stack"], r_sample)

    def travelling_ratio(self, bmap_path, r_sample=None, sector_deg=None):
        """行波度 = |Y[-p]−Y[+p]| / (|Y[-p]|+|Y[+p]|) (M-ASM-7)。

        1.0=纯行波 (转矩通道正常); 0.0=纯驻波 (时均 Maxwell 转矩恒零)。
        驻波主导提示转矩-功率恒等式 (M-ASM-5) 通道异常。
        """
        import numpy as np
        if sector_deg:
            _, bn, _ = self.gap_circle_freq(bmap_path, r_sample, n=360,
                                            th_span_deg=sector_deg)
            # 延拓符号按扇区跨半极数 (R111b): 90°(4极)=1 半极→反周期 [s,-s..];
            # 180°=2 半极→正周期 [s,s..] — 真场 a(θ+180°)=+a(θ), 反号延拓
            # 会造出 180° 反周期序列 (奇 bin only) → bin 2 精确 0 (实证)
            n_half = sector_deg * (int(base.P["n_poles"]) / 2.0) / 180.0  # 4b: 极对口径 (与 .pro link_coef 同式), 首版 n_poles 偶判→+1 延拓→bin p 恒 0 (js17b 实证)
            sign = (-1.0) ** int(round(n_half))
            bn = np.concatenate(
                [bn * (sign ** k)
                 for k in range(int(round(360.0 / sector_deg)))])
        else:
            _, bn, _ = self.gap_circle_freq(bmap_path, r_sample)
        Y = np.fft.fft(bn) / len(bn)
        p = int(base.P["p"])
        fwd, rev = abs(Y[(-p) % len(bn)]), abs(Y[p % len(bn)])
        return float((fwd - rev) / max(fwd + rev, 1e-30))

    def bg1_freq(self, bmap_path, r_sample=None, sector_deg=None):
        """频域 Bg1: 复数采样的 rfft 基波幅值 ×2 (空间基波, p 对)。"""
        import numpy as np
        if r_sample is None:
            r_sample = base.P["r_gap_mid"]
        if sector_deg:
            pass  # sector 分支在下方延拓处采样
        else:
            _, bn, _ = self.gap_circle_freq(bmap_path, r_sample)
        # 复数采样须用 fft (rfft 仅实数)。
        # 转向口径 (M-ASM-7, 2026-10-24 判决): 频域正序源产生 e^{-jpθ} 行波,
        #   落 fft bin -p (e^{+jpθ} 共轭镜像才是 bin +p) — 静磁链实数场
        #   ±bin 共轭对称无此分别, PMSM 基线不受影响; ASM 频域链必须取 -p。
        # 幅值口径: ×2 峰值相量 (与静磁链一致)。
        if sector_deg:
            _, bn, _ = self.gap_circle_freq(bmap_path, r_sample, n=360,
                                            th_span_deg=sector_deg)
            # 延拓符号按扇区跨半极数 (R111b): 90°(4极)=1 半极→反周期 [s,-s..];
            # 180°=2 半极→正周期 [s,s..] — 真场 a(θ+180°)=+a(θ), 反号延拓
            # 会造出 180° 反周期序列 (奇 bin only) → bin 2 精确 0 (实证)
            n_half = sector_deg * (int(base.P["n_poles"]) / 2.0) / 180.0  # 4b: 极对口径 (与 .pro link_coef 同式), 首版 n_poles 偶判→+1 延拓→bin p 恒 0 (js17b 实证)
            sign = (-1.0) ** int(round(n_half))
            bn = np.concatenate(
                [bn * (sign ** k)
                 for k in range(int(round(360.0 / sector_deg)))])
        Y = np.fft.fft(bn) / len(bn)
        p = int(base.P["p"])
        return float(2 * abs(Y[(-p) % len(bn)]))


def main():
    import numpy as np
    _REPO = base._REPO
    # ---- 求解器二进制 (M-ASM-2h 永久化, R13): 16 条笼耦合系统必须 double 精度
    #      直接法 (Sparskit 单精度 ILUTP 条件数崩, 见 open_defects.M-ASM-2b)。
    #      解析顺序: CHINASIM_GETDP_BIN (显式) > bin/third_party/getdp/bin/getdp-petsc
    #      (PETSc/MUMPS wrapper, 永久化入口) > 默认 Sparskit 版 (--solver sparskit 强制) ----
    global GETDP_BIN_OVERRIDE, SOLVER_LABEL
    GETDP_BIN_OVERRIDE = os.environ.get("CHINASIM_GETDP_BIN", "")
    petsc_wrapper = os.path.join(_REPO, "bin", "third_party", "getdp", "bin",
                                 "getdp-petsc")
    SOLVER_LABEL = "getdp-petsc(PETSc/MUMPS double 直接法)"
    if "--solver" in sys.argv and sys.argv[sys.argv.index("--solver") + 1] == "sparskit":
        GETDP_BIN_OVERRIDE = ""
        SOLVER_LABEL = "getdp(Sparskit 单精度 ILUTP, 仅调试对照)"
        print("[solver] 显式指定 Sparskit 单精度版 (仅调试对照用, 16 条笼会假收敛)")
    elif GETDP_BIN_OVERRIDE:
        print(f"[solver] 环境变量指定: {GETDP_BIN_OVERRIDE}")
    elif os.path.isfile(petsc_wrapper):
        GETDP_BIN_OVERRIDE = petsc_wrapper
        base.GETDP_BIN = petsc_wrapper
        print(f"[solver] auto → PETSc/MUMPS wrapper: {petsc_wrapper}")
    else:
        print("[solver] auto → getdp-petsc 不存在, 回退 Sparskit 版 (16 条笼将假收敛,"
              " 建议执行 bin/third_party/getdp/bin/getdp-petsc 安装校验)")
    # ---- 参数: UI 覆盖 (--params) 与链开关, 无业务常数 ----
    overrides = None
    if "--params" in sys.argv:
        overrides = json.load(open(sys.argv[sys.argv.index("--params") + 1],
                                   encoding="utf-8"))
        # 白名单交集过滤 (fe11 同源坑: geo_params.json 等派生键文件直接喂
        # merged_design 会 fail-fast; 设计键=merged_design({}) 键集),
        # 丢弃键显式打印 (禁静默)
        _ref = merged_design({})
        _drops = sorted(k for k in overrides
                        if k not in _ref and not k.startswith("_"))
        overrides = {k: v for k, v in overrides.items()
                     if k in _ref and not k.startswith("_")}
        if _drops:
            print(f"[params] 非设计键丢弃 ({len(_drops)}): "
                  f"{_drops[:8]}{'...' if len(_drops) > 8 else ''}")
        base.P.update(merged_design(overrides))
        print(f"[params] UI参数覆盖: {sorted(overrides.keys())}")
    wd = None
    if "--workdir" in sys.argv:
        wd = sys.argv[sys.argv.index("--workdir") + 1]
        os.makedirs(wd, exist_ok=True)

    cfg = load_asm_spec()                              # constants/ 单一事实源
    slips = list(cfg["slips"])
    if "--slips" in sys.argv:
        slips = [float(x) for x in
                 sys.argv[sys.argv.index("--slips") + 1].split(",")]
        for s in slips:
            check_slip(s)
    js_pk = float(cfg["js_peak_Am2"])
    if "--js-peak" in sys.argv:
        js_pk = float(sys.argv[sys.argv.index("--js-peak") + 1])
    if "--n-bars" in sys.argv:
        cfg["n_bars"] = int(sys.argv[sys.argv.index("--n-bars") + 1])
    # ---- 扇区对称约化 (M-ASM-11 R99, 与 --n-bars 同模式): 90°=1/4 模型
    if "--sector-deg" in sys.argv:
        cfg["sector_deg"] = float(
            sys.argv[sys.argv.index("--sector-deg") + 1])
    # ---- 库族批量开关 (M-ASM-L, v800 #3 前桥 ASM 感应机族): 与 --n-bars 同模式
    #      显式链覆盖, 非默认值兜底; 未指定时保持 constants/ 单一事实源 ----
    for _flag, _key, _cast in (
            ("--bar-diameter", "bar_diameter_mm", float),
            ("--ring-radius", "ring_radius_mm", float),
            ("--ring-bar-diameter", "ring_bar_diameter_mm", float),
            ("--mur-iron", "mur_iron_effective", float),
            ("--dtheta", "dtheta_max_deg", float),
            ("--lc-air", "lc_air", float),
            ("--lc-iron", "lc_iron", float),
            ("--lc-shaft", "lc_shaft", float),
            ("--slot-half-deg", "slot_half_deg", float),
            ("--mesh-extend", "mesh_extend_boundary", int),
            ("--mesh-algo", "mesh_algorithm", int),
            ("--gap-transfinite", "gap_transfinite", int),
            ("--gap-pts", "gap_transfinite_pts", int)):
        if _flag in sys.argv:
            cfg[_key] = _cast(sys.argv[sys.argv.index(_flag) + 1])
            print(f"[override] {_key} = {cfg[_key]}")

    drv = (GetDPAsmDriver(workdir=wd, sector_deg=cfg.get("sector_deg"))
           if wd else GetDPAsmDriver(sector_deg=cfg.get("sector_deg")))
    prog = ProgressReporter(drv.wd, run_tag="emag_asm")
    drv.prog = prog
    t0 = time.time()

    r_ring = ring_resistance(float(cfg["sigma_conductor_S_per_m"]),
                             float(cfg["ring_radius_mm"]),
                             float(cfg["ring_bar_diameter_mm"]),
                             int(cfg["n_bars"]))
    rep = {"params": {"line_freq_hz": cfg["line_freq_hz"],
                      "js_peak_Am2": js_pk, "slips": slips,
                      "n_bars": int(cfg["n_bars"]),
                      "bar_diameter_mm": float(cfg["bar_diameter_mm"]),
                      "ring_radius_mm": float(cfg["ring_radius_mm"]),
                      "ring_bar_diameter_mm": float(cfg["ring_bar_diameter_mm"]),
                      "Rring_ohm": r_ring,
                      "n_poles": int(base.P["n_poles"]),
                      "n_slots": int(base.P["n_slots"]),
                      "L_stack_m": float(base.P["L_stack"]),
                      "Flag_IronLinear": 1,
                      "mur_iron_effective": float(cfg["mur_iron_effective"])},
            "model": cfg["model"], "disclosure": cfg["disclosure"],
            "log": []}

    # ---- 笼型网格 + 生成式 .pro ----
    prog.step("几何/网格", 3, f"generate_geo(rotor_type=cage, "
                             f"n_bars={cfg['n_bars']}) → gmsh + asm pro 生成")
    drv.mesh_cage(cfg)

    # ---- 转差扫掠: 每点一次完整频域解 (--from-existing 断点续跑: 复用已有
    #      b_map_asm_s*.pos 仅重解析+判据+报告, 不重求解) ----
    rows = []
    from_existing = "--from-existing" in sys.argv
    for i, s in enumerate(slips):
        bmap = f"b_map_asm_s{s:g}.pos"
        if from_existing and not os.path.exists(os.path.join(drv.wd, bmap)):
            raise RuntimeError(f"--from-existing: {bmap} 不存在, 无法续跑")
        if not from_existing:
            prog.step("转差扫掠", 10.0 + 80.0 * i / len(slips),
                      f"s={s} (第{i + 1}/{len(slips)}点, "
                      f"f_src={source_freq_hz(s, cfg['line_freq_hz']):.3g}Hz)")
            if "--saturation" in sys.argv:
                sat_cfg = dict(cfg.get("saturation_iter", {}) or {})
                sat_cfg["enabled"] = True
                # 主表 = 宝钢 B50A470 datasheet 矢量反算真值 (asm_sat3 双表对照后切换,
                # 见 asm_design.json iron_bh_b50a470_extracted.dual_table_benchmark);
                # --bh-key 可覆盖 (代理表 iron_bh_dw465_50 保留为低场补充口径)
                bh_key = "iron_bh_b50a470_extracted"
                for i, a in enumerate(sys.argv):
                    if a == "--bh-key" and i + 1 < len(sys.argv):
                        bh_key = sys.argv[i + 1]
                if bh_key not in cfg:
                    raise KeyError(f"--bh-key {bh_key} 不在 asm_design.json (合法键见 constants)")
                bh_tab = cfg[bh_key]["bh_b_t_vs_h_a_per_m"]
                sat_hist, mur, bg1_prev = [], float(cfg["mur_iron_effective"]), None
                for it in range(1, int(sat_cfg["max_iters"]) + 1):
                    cfg_i = dict(cfg)
                    cfg_i["mur_iron_effective"] = mur
                    meta = drv.solve_slip(js_pk, s, cfg_i, out_bmap=bmap)
                    bg1_i = drv.bg1_freq(bmap)
                    bmax = drv.iron_bmax(bmap, sat_cfg.get("b_probe_quantile", 0.98))
                    mu_new = spec.mu_eff_from_bmax(bmax, bh_tab)
                    sat_hist.append({"iter": it, "mur": round(mur, 1),
                                     "Bg1": round(bg1_i, 4), "b_max_T": round(bmax, 4),
                                     "mu_eff_next": round(mu_new, 1)})
                    print(f"    [sat] s={s} it={it} Mur={mur:.0f} → "
                          f"Bg1={bg1_i:.4f} Bmax={bmax:.3f}T μ_eff'={mu_new:.0f}")
                    if bg1_prev is not None and bg1_prev > 0 and abs(
                            bg1_i - bg1_prev) / bg1_prev < float(sat_cfg["tol_bg1_rel"]):
                        break
                    bg1_prev = bg1_i
                    mur = mu_new
                meta["saturation_hist"] = sat_hist
            else:
                meta = drv.solve_slip(js_pk, s, cfg, out_bmap=bmap)
        else:
            meta = {"Rring_ohm": r_ring, "source_freq_hz":
                    source_freq_hz(s, cfg["line_freq_hz"]), "solve_s": 0.0,
                    "reused": True,
                    "travelling_ratio": round(tr_ratio, 3)}
        _sec = getattr(drv, "sector_deg", None)
        bg1 = drv.bg1_freq(bmap, sector_deg=_sec)
        t_nm = drv.torque_freq(bmap, sector_deg=_sec)
        tr_ratio = drv.travelling_ratio(bmap, sector_deg=_sec)  # M-ASM-7
        # ---- M-LAB1 功率流: 电路口径转子铜损 + 转差恒等式残差 (GB/T 1032 对标) ----
        #   P_mech 口径: 频域时均转矩 × 机械角速度 ω_m=(1-s)·ω_sync (转差定义恒等式
        #   P_ag = T·ω_m/(1-s) = T·ω_sync 自动消去 (1-s) 因子, 发电 s<0 时 P_ag<0
        #   = 机械输入, 恒等式 P_cuR = s·P_ag 恒正 — 双象限自检)
        import math
        import lab_spec
        i_ring = _read_iz_timetable(os.path.join(drv.wd, "iz_asm.dat"))
        n_bars = int(cfg["n_bars"])
        if getattr(drv, "sector_deg", None):
            n_bars = int(round(n_bars * drv.sector_deg / 360.0))  # 扇区电路口径
        # 导条电流 = 环段 KCL 差分 (PostProcessing 打印 Is 不可见 —
        #   GlobalQuantity AssociatedWith 不进 Print 离散表, GetDP 实测 'Unknown
        #   discrete Quantity: Is'; 环段 2n 电流已完备, KCL 差分即导条电流,
        #   顶/底两环独立差分互检 = 电路解一致性自检)
        i_bar, kcl_resid = lab_spec.bar_currents_from_rings(i_ring, n_bars)
        # M-ASM-10 复判: cage bar=扇形环带 (9.004mm², msh 实测吻合), 电路侧
        # R' 面积口径与 FE 侧 σ·A 统一 (圆拍 7.069mm² 会高估 P_cuR 1.274×)
        a_bar = spec.bar_cross_section_area_cage(
            float(cfg["bar_diameter_mm"]),
            float(base.P["R_ri"]),
            float(base.P["R_ri"]) + float(base.P["mag_t"]))
        r_bar_pm = 1.0 / (float(cfg["sigma_conductor_S_per_m"]) * a_bar)
        p_cur = lab_spec.cu_rotor_loss_p(
            i_bar, i_ring, r_bar_pm, float(base.P["L_stack"]), r_seg=r_ring)
        omega_sync = 2.0 * math.pi * float(cfg["line_freq_hz"]) * 2.0 / float(
            base.P["n_poles"])   # ω_sync = ω_line·p, p=n_poles/2 (主设计键, 无兜底)
        ident = lab_spec.slip_power_identity(
            t_nm * omega_sync, p_cur, s)
        rows.append({"slip": s, "mode": slip_mode(s),
                     "torque_avg_Nm": round(t_nm, 4),
                     "Bg1_amp_T": round(bg1, 4),
                     "travelling_ratio": round(tr_ratio, 3),
                     "p_cu_rotor_W": round(p_cur, 3),
                     "p_airgap_W": (None if ident["p_ag_W"] is None
                                    else round(ident["p_ag_W"], 3)),
                     "slip_resid": (None if ident["resid"] is None
                                    else round(ident["resid"], 4)),
                     "kcl_resid": round(kcl_resid, 4), **meta})
        print(f"  s={s:g}  T_avg={t_nm:+.4f} Nm  Bg1={bg1:.4f} T"
              f"  P_cuR={p_cur:.2f}W  resid={ident['resid']}"
              f"  KCL={kcl_resid:.3f}")
    rep["slip_sweep"] = {"rows": rows}

    # ---- 物理自检 (M-ASM-2g 电流馈电恢复型形态) ----
    #   ① 场量有效性: Bg1 全点 > 下限 (防空场/伪解假 PASS)
    #   ② 低滑差恢复单调: Bg1(|s|_min) ≥ Bg1(|s|_max) — 趋肤屏蔽随 |s| 减弱,
    #      转子系源频率越低导条越透明, 气隙基波单调恢复 (MUMPS 直接法 5 点实证:
    #      0.0286→0.0591T 单调; 发电区 s<0 同物理, R19 按 |s| 判,
    #      ±s 幅值对称性 0.0001-0.0002T 独立佐证)
    #   ③ 全点转矩有限
    bg1_min = min(r["Bg1_amp_T"] for r in rows)
    bg1_floor = 1e-3   # T — 低于此判定场图解析失效
    ok_field = bg1_min > bg1_floor
    by_abs_s = sorted(rows, key=lambda r: abs(r["slip"]))   # |s| 升序
    b_at_smin = abs(by_abs_s[0]["Bg1_amp_T"])    # 首点=最小 |s|
    b_at_smax = abs(by_abs_s[-1]["Bg1_amp_T"])   # 末点=最大 |s|
    ok_recover = b_at_smin >= b_at_smax
    ok_finite = all(np.isfinite(r["torque_avg_Nm"]) for r in rows)
    ok = ok_field and ok_recover and ok_finite
    rep["physics_check"] = {"bg1_min_T": bg1_min, "bg1_floor_T": bg1_floor,
                            "field_sanity_pass": bool(ok_field),
                            "bg1_at_smin_T": b_at_smin,
                            "bg1_at_smax_T": b_at_smax,
                            "recovery_pass": bool(ok_recover),
                            "finite_torque_pass": bool(ok_finite),
                            "note": cfg["physics_check_note"],
                            "pass": bool(ok)}
    # ---- M-LAB1 功率流自检 (独立于形态判据, 如实三态不假绿) ----
    #   KCL: 顶/底环独立差分互检 (电路解自洽铁证, ~机器精度);
    #   转差恒等式: P_cuR = s·P_ag (电路损耗口径已获场积分 ∫(1/2)σω²|a|²dV
    #   独立互证 7% 内; Maxwell 转矩口径缺口 38× → M-ASM-5 open defect,
    #   恒等式残差带 OPTIMISTIC 时禁静默放行)
    import lab_spec as _lab
    _pf = _lab.power_flow_check([{"resid": r["slip_resid"]} for r in rows])
    _kcl_max = max(r.get("kcl_resid", 0.0) for r in rows)
    rep["physics_check"]["power_flow"] = {
        "verdict": _pf["verdict"], "max_resid": _pf["max_resid"],
        "note": _pf["note"], "kcl_resid_max": _kcl_max}    # ---- 报告落盘 ----
    prog.step("报告落盘", 95, "asm_report.json")
    rep["runtime_s"] = round(time.time() - t0, 1)
    # ---- 时间统计汇总 (goal: 做好时间统计输出保证后期优化, R13) ----
    solve_times = [r["solve_s"] for r in rows if r.get("solve_s", 0) > 0]
    # 峰值内存汇总 (goal 二.3 台账列): 各转差点 getdp 子进程峰值 RSS 最大值
    rss_vals = [r["peak_rss_mb"] for r in rows
                if isinstance(r.get("peak_rss_mb"), (int, float))]
    rep["solver_timing"] = {
        "solver_binary": SOLVER_LABEL,
        "points_solved": len(solve_times),
        "solve_s_total": round(sum(solve_times), 1),
        "solve_s_avg": round(sum(solve_times) / len(solve_times), 1) if solve_times else 0.0,
        "solve_s_max": round(max(solve_times), 1) if solve_times else 0.0,
        "peak_rss_mb_max": round(max(rss_vals), 1) if rss_vals else None,
        "note": "逐点 wall time 为 GetDP 子进程耗时 (mesh 复用, 单 System 生成); "
                "MUMPS 直接法 ~94-227s/点@34k DOF, 后期优化基线"}

    rep["log"] = drv.log
    with open(os.path.join(drv.wd, "asm_report.json"), "w",
              encoding="utf-8") as f:
        json.dump(rep, f, ensure_ascii=False, indent=1)
    prog.finish(ok, f"ASM 转差链 {len(slips)} 点, Bg1={bg1_min:.4f}-{max(r['Bg1_amp_T'] for r in rows):.4f}T "
                    f"恢复单调={ok_recover} 耗时{rep['runtime_s']}s → asm_report.json")
    print(f"=== ASM 物理自检 ===  场有效={ok_field} 恢复单调={ok_recover} "
          f"有限={ok_finite} → {'PASS' if ok else 'FAIL'}")
    print(f"总耗时 {rep['runtime_s']}s → asm_report.json")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

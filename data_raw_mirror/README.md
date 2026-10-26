# data/raw — 数据地图（论文数字 ↔ 归档工件溯源表）

每个正文数字的来源文件。铁律：正文数字 100% 可溯源，无手输。

## Case A — 求解器预算中止（§Case A, Fig.1, tab:cases 行 A）

| 工件 | 内容 | 支撑的正文数字 |
|---|---|---|
| `case1_cross_solver/fe17_chi2_mumps_fail_run.out` | **χ=2.0 bar 畸形几何 MUMPS 一手 stdout**（37.0s 全 nan FAIL, 自检 False→FAIL） | `-19` 失败形态、all-NaN 传导 |
| `case1_cross_solver/front2_fe17_params.json` / `front2_fe18_params.json` | χ=2.0 (bar 12 畸形) 与 χ=3.067 (bar 9.2 归一) 参数快照 | 复现参数 |
| ⏳ PENDING | **串行 SuperLU 同矩阵一手日志**（历史 campaign fe19/fe20 真跑: 26s/残差 2.28e-5/T=0.7477/Bg1=0.0158, 一手 stdout 未及时归档已随 .tmp 清理丢失; 重跑脚本 `../src/run_e7c_true_cross.sh` 待低负载窗口执行后补档） | 摘要交叉数字的最终一手溯源 |
| `e7_protocol_fail/step1_mumps.log` / `step4_superlu.log` | ⚠️ 更正: 实为 χ=2.0 bar 归一 **PASS 演练** 双跑 (T=0.2495/Bg1=0.0363, 均 MUMPS 通道; 归一后健康解用于 A4 披露语义验证), **非交叉求解器对** | E7b 归一 PASS 事实 |
| `e7_protocol_fail/E7B_FINDING.md` | bar 归一精化：χ=2.0 归一 PASS → 失败带 (2.0, 3.067] | Fig.1 失败带位置 |
| `e7_protocol_fail/fe17_params_x2.json` | χ=2.0 参数快照（白名单键） | 复现参数 |
| `mumps_infog19_semantics.md` | **MUMPS 5.8.2 源码级证据**（dfac_b.F:283：-19=ICNTL(23) 预算超限，非奇异码） | M1 修正后的标签定性 |
| `case1_nocir/nocir_bisection_from_evidence.md` | 纯场二分（电路全注释仍失败）→ 排除电路耦合 | 排除矩阵 (i) |
| `case1_scale_ladder/asm_full16d_report.json` | χ=1 基线 11 点结构化报告（solver_timing: 62.6 s/223 MB） | 基线 `0.0579`(s=1)/`0.1141`(s=0.05) |
| `e7_protocol_pass/asm_report.json` + `step1_mumps_pass.log` | PASS 分支演练（MUMPS 首解即成功） | 协议 PASS 分支 |
| `case_evidence_from_asm_design.json` | M-ASM-15/16/17/18 结构化 evidence（campaign 期归档） | Case B χ=1.5/1.75 数字、M-ASM-18 物理定性 |

## Case B — 中间尺度半场（§Case B, Fig.1 中 χ=1.5 点）

χ=1.5/1.75 数字（0.0663/0.1141/0.1161，同工况 s=0.05, js=3 A/mm²）来自
归档实验 fe16/fe15/fe14b，其结构化证据固化在
`case_evidence_from_asm_design.json` 的 M-ASM-15 字段（阈值谱 evidence）。

## Case C — 扇区均匀符号违例（§Case C, Fig.3, tab:cases 行 C）

| 工件 | 内容 | 支撑的正文数字 |
|---|---|---|
| `case2_sector/case2_numbers_from_opslog.md` | js17b 扇区 vs asm_full16d 全模型 | `T 0.0605/0.0393 = ×1.54`、`Bg1 0.0833/0.0579 = ×1.44`、KCL=0.000 |
| `case2_sector/…` 同目录 | 阶次相容性分析（12k 齿谐波 vs 90°） | uniform-sign condition |

## 重现脚本

`../src/run_e7_protocol.sh`（FAIL 分支五步）/ `../src/run_e7a_pass_branch.sh`
（PASS 分支）；审计参考实现 `../src/audit_protocol.py`（doctest 驱动，
含两个合成夹具：NaN 传导与交叉失配）。

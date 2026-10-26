# Case 1 纯场二分 — 真数字 (来源: M-ASM-15 evidence, motor_mag_asm_nocir.pro)

- 手法: GlobalTerm 70 行注释 + GlobalEquation 括号平衡整块注释 (If(0) 包裹法语法碎教训)
- 结果: ×2.0 与 ×3.067 纯场 (无电路自由度) 均 MUMPS -19
- 结论历史: 当时判"奇异在场块 (a 逐点+ur 区域标量)" → 后被 SuperLU 交叉解 (E2)
  推翻为 "MUMPS 求解器特定失败" — 二分实验保留为排除链一环
- 工程教训: 手改 .pro 单行注释安全; GlobalEquation 整块需括号平衡行级注释

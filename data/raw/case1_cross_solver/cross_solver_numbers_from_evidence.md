# Case 1 交叉求解 — 真数字 (来源: constants/asm_design.json M-ASM-15 evidence 字段, git R136-R140)

## 同矩阵双求解器对照 (fe14 几何 R_si 92mm ×3.067, 16 条笼, js=3e6, s=1, 全模型 116476 节点/232078 三角)
| 项 | MUMPS (PETSc wrapper) | SuperLU (serial, PETSc --download-superlu) |
|---|---|---|
| 求解结果 | ERROR -19 (singular, 假奇异) | 收敛, residual 2.28e-5 |
| 解的量值 | 无解 | T=0.7477 Nm, Bg1=0.0158 T, KCL=0.000 |
| 求解耗时 | 失败前 ~分钟级 (OOC ICNTL22=1+23=8000) | 26 s |
| 全流程 | — | 125.9 s PASS |

## 排除矩阵摘要 (E4, 每行=一个假设的否证)
- OOC 越内存墙: -13(alloc) → -19, 排除"内存不足假象"
- ICNTL(6)=77+(8)=77 缩放: 无效, 排除"数值缩放病"
- ICNTL(24)=1+CNTL3 null pivot 检测 / ICNTL(11)=6 / ICNTL(13)=6 强制 stdout: 全无输出 (wrapper 屏蔽), 定位通道缺失
- 纯场二分 (GlobalTerm 70 行+GlobalEquation 整块注释): ×2/×3.067 仍 -19 → 奇异在场块非电路 → 后被 E2 推翻为"假奇异"
- 条数无关: 16/28 条均崩; bar-ratio 归一 (fe18 bar 9.2@R_si 92): 仍 FAIL → 排除几何比例混杂
- 网格质量: min 面积 2.4e-8 无退化单元; 26 组全非空 → 排除网格/组
- 阈值谱: ×1.0 PASS (Bg1=0.1141) / ×1.5 病态半场 (0.0663=58%) / ×1.75 PASS (0.1161, 自洽 1.7%) / ×2.0 崩 (37s 全 nan) / ×3.067 崩 → 崩塌 ∈(1.75,2.0]

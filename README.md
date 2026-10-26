# Paper 11 — Direct-Solver Silent Failure in Magnetoquasistatic FEM:
# False Singularity Diagnosis and a Cross-Solver Verification Protocol

- 状态: S1 选题价值门禁进行中 (2026-10-26 启动)
- 选题: P8-1 主攻 (磁-电路耦合 FEM 直接法静默失败/假奇异实证 + MUMPS vs SuperLU
  交叉验证协议 + 尺度阈值谱) + P8-2 并入 (周期对称约束齿谐波污染定量机理)
- 数据源: ChinaSimStdio R136-R219 五轮悬案全链真数据 (真数据铁律天然满足)
- 求解器口径: GetDP 4.0.0 + PETSc 3.25.5 (3.5.0 源码树已删除, 一律 4.0.0)
- 自包含铁律: 代码/数据/脚本/结果全在本目录; 三方库只读引用 bin/third_party

## 主稿口径 (2026-10-26 R3, 路线 A 重定位)
- 主稿: docs/manuscript/main_compel.tex (elsarticle 单栏, COMPEL 目标) —
  "Application-Level Silent Failures in Electrical-Machine FE: Taxonomy and Audit Protocol"
- 编译: cd docs/manuscript && pdflatex main_compel.tex (×2) → 7 页 0 error 4 图
- IEEE 双栏旧版: main.tex (保留存档)
- 数据唯一来源: data/raw/ (六目录), 图源: figure/fig1-4 (src/make_figures.py 生成)
- 查新: 00_立项卡.md (三轮: S1 价值门禁 34 篇/四结论查新/SDC 与 Result-Verification 邻居差异化)

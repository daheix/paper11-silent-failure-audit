#!/usr/bin/env bash
# M12 — Case B (chi=1.5) bar 归一对照复跑
# 目的: 判定 0.0663 (58% 半场) 是尺度条件数效应还是 bar 畸形混杂
#   fe16 口径: R_si 30 未放大 + bar 9.0 (ratio 0.30 畸形) -> Bg1=0.0663
#   本脚本:   R_si 45 (真 chi=1.5) + bar 4.5 (ratio 0.10 归一) -> Bg1=?
# 判读: 归一后 ~0.17 (×1.5 of 0.1141) = 尺度效应成立(F2 keep, fe16 归混杂)
#        归一后仍 ~0.066 = F2 尺度病成立; 介于其间 = 立项再议
# 前置: load < 3.5
set -euo pipefail
GETDP_DIR=/home/wsc/wsc/ChinaSimStdio/src/plugins/motor_workbench/simulation/getdp
OUT="$(cd "$(dirname "$0")/.." && pwd)/data/raw/case1_scale_ladder/m12_chi15_normalized"
LOAD=$(awk '{print $1}' /proc/loadavg)
OK=$(python3 -c "print(1 if float('$LOAD')<3.5 else 0)")
if [ "$OK" != "1" ]; then echo "load=$LOAD >= 3.5, 窗口未开, 跳过"; exit 0; fi
mkdir -p "$OUT"; cd "$GETDP_DIR"
WD=../../.tmp/m12_chi15_norm
rm -rf "$WD"; mkdir -p "$WD"
python3 run_emag_asm.py --workdir "$WD" --slips 0.05 \
  --js-peak 3000000 --n-bars 16 --bar-diameter 4.5 --ring-radius 52.5 \
  --ring-bar-diameter 6.0 --dtheta 2.5 --solver mumps \
  2>&1 | tee "$OUT/run.log"
cp "$WD/asm_report.json" "$OUT/"
echo "M12 完成: $OUT"

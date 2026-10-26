#!/usr/bin/env bash
# E7a 协议 PASS 分支演练 - js18 口径小几何 (全模型 16 条, MUMPS 首解成功路径)
set -euo pipefail
GETDP_DIR=/home/wsc/wsc/ChinaSimStdio/src/plugins/motor_workbench/simulation/getdp
BASE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$BASE/data/raw/e7_protocol_pass"; mkdir -p "$OUT"
cd "$GETDP_DIR"
unset PETSC_OPTIONS
python3 run_emag_asm.py --workdir ../../.tmp/e7a_pass --slips 1.0 --solver mumps \
  --js-peak 3000000 --n-bars 16 --bar-diameter 3.0 --ring-radius 35.0 \
  --ring-bar-diameter 4.0 --dtheta 2.5 2>&1 | tee "$OUT/step1_mumps_pass.log"
cp ../../.tmp/e7a_pass/asm_report.json "$OUT/"
echo "E7a done: $OUT"

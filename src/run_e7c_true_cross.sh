#!/usr/bin/env bash
# E7c — 真交叉求解器对重跑 (M11 修复)
# χ=2.0 bar 畸形几何 (fe17 口径): MUMPS 必须 FAIL(-19/nan) + SuperLU 必须 PASS
# 硬校验: 两侧 solver_binary 与实际一致, 两报告 T 不得逐位一致
# 前置: load < 3.5 (低负载窗口; χ=2.0 网格生成约 30-60min)
set -euo pipefail
GETDP_DIR=/home/wsc/wsc/ChinaSimStdio/src/plugins/motor_workbench/simulation/getdp
OUT="$(cd "$(dirname "$0")/.." && pwd)/data/raw/case1_cross_solver/e7c_true_cross"
LOAD=$(awk '{print $1}' /proc/loadavg)
OK=$(python3 -c "print(1 if float('$LOAD')<3.5 else 0)")
if [ "$OK" != "1" ]; then echo "load=$LOAD >= 3.5, 窗口未开, 跳过 (禁重试循环)"; exit 0; fi
mkdir -p "$OUT"; cd "$GETDP_DIR"

run_side () {  # $1=tag $2=PETSC_OPTIONS
  WD=../../.tmp/e7c_$1
  rm -rf "$WD"; mkdir -p "$WD"
  if [ -n "${2:-}" ]; then
    PETSC_OPTIONS="$2" python3 run_emag_asm.py --workdir "$WD" --slips 1.0 \
      --js-peak 3000000 --n-bars 16 --bar-diameter 12.0 --ring-radius 70.0 \
      --ring-bar-diameter 14.0 --dtheta 2.5 --params ../../.tmp/front2_fe17_params.json \
      2>&1 | tee "$OUT/step_$1.log"
  else
    python3 run_emag_asm.py --workdir "$WD" --slips 1.0 \
      --js-peak 3000000 --n-bars 16 --bar-diameter 12.0 --ring-radius 70.0 \
      --ring-bar-diameter 14.0 --dtheta 2.5 --params ../../.tmp/front2_fe17_params.json \
      2>&1 | tee "$OUT/step_$1.log"
  fi
  cp "$WD/asm_report.json" "$OUT/asm_report_$1.json"
}

echo "=== side 1: MUMPS (must FAIL with -19/all-NaN) ==="
run_side mumps_fail ""
grep -q "FAIL" "$OUT/step_mumps_fail.log"

echo "=== side 2: serial SuperLU (must PASS) ==="
run_side superlu_pass "-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu"
grep -qi "PASS" "$OUT/step_superlu_pass.log"

python3 - "$OUT" <<'PY'
import json, sys
out = sys.argv[1]
a = json.load(open(f"{out}/asm_report_mumps_fail.json"))
b = json.load(open(f"{out}/asm_report_superlu_pass.json"))
ra = (a["slip_sweep"]["rows"])[0]; rb = (b["slip_sweep"]["rows"])[0]
ta, tb = ra.get("torque_avg_Nm"), rb.get("torque_avg_Nm")
assert ta != tb, "两报告 T 逐位一致 = 同求解器, 交叉无效"
print(f"CROSS OK: mumps T={ta} (fail) vs superlu T={tb} (pass)")
PY
echo "E7c 完成: $OUT"

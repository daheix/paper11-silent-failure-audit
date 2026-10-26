#!/usr/bin/env bash
# E7b 协议 FAIL 分支 - fe17 口径 (x2.0 尺度: R_si 60, M-ASM-15 崩塌带内)
# 共存模式: nice -n 19 不夺 paper08 频扫算力; 门禁放宽为 load<7 (nice 自缚)
set -uo pipefail
GETDP_DIR=/home/wsc/wsc/ChinaSimStdio/src/plugins/motor_workbench/simulation/getdp
BASE="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$BASE/data/raw/e7_protocol_fail"; mkdir -p "$OUT"
cd "$GETDP_DIR"
load1=$(awk '{print $1}' /proc/loadavg)
python3 - "$load1" <<'PY' || { echo "load>=7, 禁止启动"; exit 2; }
import sys; sys.exit(0 if float(sys.argv[1]) < 7 else 1)
PY
# x2.0 params (asm_full16d 基线 x2, 白名单键口径, R137 阈值谱 fe17 复刻)
python3 - "$OUT" <<'PY'
import json, sys, os
base = json.load(open('../../.tmp/asm_full16d/geo_params.json'))
s = 2.0
for k in ('R_si','R_rm','R_ri','R_so','R_shaft','slot_depth'):
    base[k] = round(base[k]*s, 6)
base['mag_t'] = round(base.get('mag_t', 0.004)*s, 6) if base.get('mag_t') else base.get('mag_t')
out = sys.argv[1]
json.dump(base, open(os.path.join(out,'fe17_params_x2.json'),'w'), indent=1)
print('params x2:', {k: base[k] for k in ('R_si','R_rm','R_so','n_slots','n_poles')})
PY
# 协议 step1: MUMPS 首解 (期望失败: nan / -19 系)
unset PETSC_OPTIONS
python3 run_emag_asm.py --workdir ../../.tmp/e7b_fail --slips 1.0 --solver mumps \
  --js-peak 3000000 --n-bars 16 --bar-diameter 6.0 --ring-radius 70.0 \
  --ring-bar-diameter 8.0 --dtheta 2.5 --params "$OUT/fe17_params_x2.json" \
  2>&1 | tee "$OUT/step1_mumps.log"
cp ../../.tmp/e7b_fail/asm_report.json "$OUT/" 2>/dev/null || true
# 协议 step4: SuperLU 交叉解
export PETSC_OPTIONS="-ksp_type preonly -pc_type lu -pc_factor_mat_solver_type superlu"
python3 run_emag_asm.py --workdir ../../.tmp/e7b_fail_sl --slips 1.0 --solver mumps \
  --js-peak 3000000 --n-bars 16 --bar-diameter 6.0 --ring-radius 70.0 \
  --ring-bar-diameter 8.0 --dtheta 2.5 --params "$OUT/fe17_params_x2.json" \
  2>&1 | tee "$OUT/step4_superlu.log"
cp ../../.tmp/e7b_fail_sl/asm_report.json "$OUT/asm_report_superlu.json" 2>/dev/null || true
echo "E7b done: $OUT"

#!/usr/bin/env bash
#SBATCH --job-name=dgs-fixed
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --array=0-59%40
#SBATCH --output=results/logs/simulation_%A_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"${PYTHON_BIN:-python}" -u run_study.py --task-index "$SLURM_ARRAY_TASK_ID" --slim "${SLIM_BIN:-slim}"

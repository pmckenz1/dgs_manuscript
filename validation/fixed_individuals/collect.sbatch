#!/usr/bin/env bash
#SBATCH --job-name=dgs-fixed-figure
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=results/logs/figure_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"${PYTHON_BIN:-python}" -u run_study.py --collect-only
"${PYTHON_BIN:-python}" -m jupyter nbconvert --to notebook --execute --inplace simulation_study_consolidated_results.ipynb

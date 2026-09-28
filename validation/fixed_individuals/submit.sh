#!/usr/bin/env bash
# Pass site-specific sbatch options, e.g. bash submit.sh --account=LAB --partition=QUEUE
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
mkdir -p results/logs
simulation=$(sbatch --parsable "$@" simulations.sbatch)
simulation=${simulation%%;*}
printf 'simulation\t%s\n' "$simulation" >> results/slurm_jobs.tsv
figure=$(sbatch --parsable "$@" --dependency="afterok:$simulation" collect.sbatch)
figure=${figure%%;*}
printf 'figure\t%s\n' "$figure" >> results/slurm_jobs.tsv
printf 'Simulation array: %s\nDependent figure job: %s\n' "$simulation" "$figure"

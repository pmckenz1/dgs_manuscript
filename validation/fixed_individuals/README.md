# Fixed-individual likelihood heatmap

This directory contains the full workflow: simulate populations in SLiM, sample fixed individuals, calculate DGS likelihoods, and plot the saved results. Install `selfdgs` separately.

These single-realization analyses use fixed individuals and a shared multilocus
pedigree. The manuscript simulation check uses independent populations and
freshly sampled individuals at each locus. Nested sample sizes are
paired observations, not independent population replicates. The larger-sample
heatmap also changes sample size and uses a display floor of −200, compared
with −80 in the preprint; use a common scale for direct visual comparisons.
At N=100, n=100 samples the entire population. None of these likelihood widths
is a calibrated uncertainty interval.

## Plot the supplied results

Open `simulation_study_consolidated_results.ipynb` and run its two code cells. The figure has sample-size rows n=20, 50, 100 and census-size columns N=100, 250, 500. `CLIP_FLOOR = -200.0` controls the shared color scale.

The notebook reads `results/likelihoods.csv` and writes `results/fixed_individuals_likelihood_heatmap.png` and `.pdf`. Plotting requires only Python, NumPy, pandas, Matplotlib, and Jupyter; it does not require SLiM or selfdgs.

Alternatively, execute it from this directory without opening a notebook editor:

```bash
python -m jupyter nbconvert --to notebook --execute --inplace simulation_study_consolidated_results.ipynb
```

## Reproduce from simulations

Use Python 3.12, SLiM **5.2** on PATH, and an installed `selfdgs` package (version **0.1.0** was used). `requirements.txt` records the supporting Python package versions used for the supplied results:

```bash
python -m pip install -r requirements.txt
python -m pip install git+https://github.com/pmckenz1/selfdgs.git
python -c "import selfdgs; print(selfdgs.__version__)"
slim -v
```

From this directory, run all 60 populations sequentially:

```bash
python run_study.py
```

Then run the plotting notebook. `--slim /path/to/slim` selects a binary outside PATH. Compatible completed simulations are reused on subsequent runs. Each population is fitted at all three sample sizes. Once all 60 finish, their likelihood tables replace `results/likelihoods.csv`. To reproduce into a separate directory without replacing the supplied table:

```bash
python run_study.py --output results/reproduced
```

The new table is `results/reproduced/likelihoods.csv`; set `RESULTS = HERE / "results/reproduced"` in the notebook to plot it.

Other useful commands:

```bash
python run_study.py --describe                # all task parameters and seeds
python run_study.py --task-index 0            # one population, all three samples
python run_study.py --analyze-only            # refit existing population VCFs
python run_study.py --collect-only            # combine all 60 completed fit tables
```

Tasks 0–19 use N=100, 20–39 use N=250, and 40–59 use N=500; within each block selfing increases from 0 to 0.95. Do not run overlapping tasks writing to the same output directory.

## Optional Slurm execution

Activate the Python environment with selfdgs and ensure SLiM 5.2 is on PATH. Then:

```bash
bash submit.sh --account=YOUR_ACCOUNT --partition=YOUR_PARTITION
```

The array runs at most 40 tasks simultaneously, each requesting one CPU, 8 GB, and four hours. A dependent job collects the likelihoods and executes the plotting notebook only after every array task succeeds. `PYTHON_BIN` and `SLIM_BIN` can override executable names. Logs and job IDs are written under `results/`. Adjust the site-specific account, partition, and resource requests for your cluster.

## Methods and files

- `fixed_individuals.slim`: the actual multilocus simulation model. All loci share a population pedigree. One nonmutating spacer with recombination probability 0.5 separates adjacent loci. The full final population is exported in a fixed VCF-column order.
- `run_study.py`: seed generation, simulation execution, VCF parsing, nested individual selection, DGS counts, likelihood calculation, and collection.
- `simulation_study_consolidated_results.ipynb`: the short plotting notebook.
- `simulations.sbatch`, `collect.sbatch`, `submit.sh`: optional cluster wrappers; the sequential runner works without Slurm.
- `results/likelihoods.csv`: 3,600 candidate evaluations, covering 180 fits from 60 populations.
- `results/fixed_individuals_likelihood_heatmap.png` and `.pdf`: exported notebook figures.

The CSV has one row per census size, true selfing rate, sample size, and candidate
selfing rate. `census_size` and `n_diploids` count diploid individuals; `true_s`
and `candidate_s` are selfing probabilities. `replicate` is always 0 and `n_loci`
is always 10,000. `loglik` is the unnormalized summed log likelihood;
`delta_loglik` subtracts the maximum over candidates within each fit. Values in
the CSV are not clipped. Compare relative likelihoods within a fit, not absolute
log likelihoods between sample sizes or populations.

Raw population VCFs and intermediate files are not bundled. A simulation rerun
creates `results/populations/N{N}/s_{rate}/` containing `population.vcf`,
`histories.csv`, `simulation.json`, SLiM logs, `sampled_individuals.csv`,
`observed_spectra.csv`, `analysis.json`, and each population's `likelihoods.csv`.
These retain genotypes, sample membership, observed counts, seeds, versions,
and code hashes for inspecting the regenerated analysis. `--analyze-only` and
`--collect-only` require those regenerated files; the supplied aggregate CSV
alone is sufficient for plotting.

Each population contains 10,000 loci of 1,000 bp; mutation is 1e-7 per base per generation, within-locus recombination is 5e-8, and burn-in is 10 × census size generations. True and candidate selfing rates are 0–0.95 in increments of 0.05. Census sizes are 100, 250, 500; there is one population realization per rate/census combination. The same sampled individuals are held fixed across all loci. Samples use the first 20, 50, or 100 entries of a random ordering of 100 distinct individuals.

For every segregating biallelic site, the DGS counts individuals carrying 0, 1, or 2 alternate alleles. Both monomorphic configurations are excluded. The likelihood is the sum of observed configuration counts times log analytical selfdgs probabilities, conditional on segregating sites, unfolded with reference polarization. The model probabilities are cached across populations. The notebook normalizes each dataset's log likelihood to its maximum and clips only for display.

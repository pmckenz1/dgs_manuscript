# Manuscript simulation check

This directory contains the independent-population simulation analyses used
for Figures 2–3. The 20 rate notebooks cover selfing rates
0.00–0.95 in steps of 0.05, census sizes N=100, 250, and 500, and samples of
n=3 and n=6 diploid individuals.

Each parameterization pools 10,000 independently simulated loci/populations. At each
locus, n=3 and n=6 are separate draws from the same final population. Individuals
are resampled at the next locus. There is one fitted dataset per condition.
These analyses check marginal DGS expectations. Performance and uncertainty
for individuals held fixed across a genome require a different sampling design;
see the [fixed-individual analyses](../fixed_individuals/README.md).

## Use the saved results

Start Jupyter in this directory. With numpy, pandas, matplotlib, and jupyter
installed, run `simulation_study_consolidated_results.ipynb` to combine the
120 saved sample-specific fits and plot the recovery summaries and likelihood
heatmap. It reads CSV files under `results/joint_s*/N*/n*/` and writes tables
and PNG/PDF figures to `results/consolidated/`. It does not require SLiM or
`selfdgs`.

The individual `simulation_study_s0-*.ipynb` notebooks also require an installed
`selfdgs` package.
Leave `EXECUTE_SLIM=False` to load saved results and redraw the per-rate figures.
The saved combined VCFs are included because the rate notebooks check for them
when loading outputs, even though plotting uses the saved tables.

## Regenerate simulations

Install SLiM 5.2 and `selfdgs`.
Set `EXECUTE_SLIM=True` in the desired rate notebook and run its cells from
this directory. This regenerates the simulations and replaces that rate's
saved results. Run the consolidated notebook after all desired rates finish.

The notebooks specify mutation rate 1e-7, within-locus recombination rate 5e-8,
1,000 bp per locus, and burn-in of 10 times census size. Parameters are saved
in each rate's `simulation_parameters.json`, and `N*/locus_seeds.csv` records
the seeds shared by the two sample-size draws. Each `n3/` or `n6/` directory
contains observed DGS counts, likelihoods, fit metadata, a validation summary,
and a combined VCF. Saved tables, VCFs, and exported figures are included.

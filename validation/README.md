# Manuscript validation analyses

## Validating theory with sims

[`manuscript_simulation_check/`](manuscript_simulation_check/README.md) contains the analyses used in
Figures 2–3. Each fitted dataset pools 10,000 independently simulated
locus/population realizations with fresh individual samples at each locus.
There is one fitted dataset per parameterization.
The goal here is to validate marginal DGS expectations, not to estimate performance or
uncertainty for individuals held fixed across the genome, see below.

## Fixed-individual analyses

See [`fixed_individuals/README.md`](fixed_individuals/README.md) for the
fixed-individual likelihood heatmaps and complete reproduction instructions. 
These simulations retain the sampled individuals' shared selfing histories
across loci and fit the marginal DGS likelihood.

The study compares nested samples of 20, 50, and 100 diploid individuals at
census sizes 100, 250, and 500 across 20 selfing rates (0–0.95 by 0.05).
Each population has 10,000 loci sharing a pedigree, and the same sampled
individuals are used across loci. There is one population per rate/census
combination; the panels do not summarize variability across replicates.

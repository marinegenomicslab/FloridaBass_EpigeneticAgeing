# Weber et al. 2025 (Submitted to Ecology and Evolution)

Data and scripts/code associated with "Nonlethal, epigenetic age estimation in a freshwater sportfish, Florida bass (*Micropterus salmoides*)" D.N. Weber et al. 2025. Submitted to Ecology and Evolution.

## Associated data can be found in the "Data" folder, including:
- "RADmeth_treated_filtered.zip" contains methylated and total read counts for all individuals and all CpG sites post-filtering for depth and coverage.
- "Metadata.csv" contains the metadata associated with all individuals included in the present study (including sample IDs and otolith-derived ages).

## Associated code can be found in the "Scripts" folder, including:

- "BAYES_GLM.r" contains code associated with the Bayesian GLM run with 4,000 warm-up/sampling iterations.
- "95% CI's.Rmd" contains the code necessary to calculate 95% confidence intervals.
- "Elastic Net Loop.r" is the custom R script used to randomly sample loci and run X iterations of elastic net.

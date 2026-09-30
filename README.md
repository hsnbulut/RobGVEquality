# Minimal Dataset and Reproducibility Files

This archive supports the results reported in the manuscript
“A Robust Permutation Test for Equality of Generalized Variances Based on MCD
Scatter Estimation”.

## Related manuscript

Bulut, H.; Gürlen, Ü. A Robust Permutation Test for Equality of Generalized
Variances Based on MCD Scatter Estimation. Unpublished manuscript, 2026.

## Contents

- `data/final_R2000_B1999/`: one RDS file for each of the 50 simulation
  scenarios. The files contain the simulation outputs and scenario-level
  results from the final run, using 2,000 Monte Carlo replications and 1,999
  permutation replicates.
- `data/realdata_crabs_results.csv`: results for the original and controlled
  rowwise-contaminated crab morphology data.
- `code/`: R code used to generate the data, run the competing procedures,
  summarize simulation results, and analyze the crab data.

## Data description

The simulation study considers independent multivariate samples under clean
Gaussian models, clean alternatives, and directional rowwise contamination.
The scenarios vary the number of groups, dimension, sample size, generalized-
variance signal, and contamination fraction. The response columns in the
real-data file report the test statistic, p-value, decision at alpha = 0.05,
and the number of successful and failed MCD permutations.

The crab morphology data are not redistributed in this archive because they
are third-party data. They are publicly available through the `crabs` dataset
in the R package `MASS` and can be loaded with:

```r
data(crabs, package = "MASS")
```

## Software

The analysis requires R and the packages `MASS`, `rrcov`, and `MVTests`.
The proposed procedure is implemented as `RobPer_GVTest()` in `MVTests`.

The RDS files can be read with `readRDS()`. The CSV file can be read with
`read.csv()`.

## Applying the proposed method to real data

The proposed test can be applied to a new multivariate data set using the
`RobPer_GVTest()` function in the `MVTests` R package:

```r
library(MVTests)

fit <- RobPer_GVTest(
  x = X,
  group = group,
  B = 1999,
  alpha = 0.75,
  seed = 123
)

fit$statistic
fit$p.value
summary(fit)
```

The package is available at:
https://github.com/hsnbulut/MVTests

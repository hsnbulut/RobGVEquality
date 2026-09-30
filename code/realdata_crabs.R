#!/usr/bin/env Rscript

## Real-data illustration for the generalized-variance tests.
## Primary analysis: BM, BF, and OM groups from MASS::crabs.
## The contaminated analysis is a pre-specified rowwise perturbation,
## not a data-driven search for a desired result.

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1L) args[[1L]] else "realdata_crabs_results.csv"
B <- as.integer(Sys.getenv("B_REAL", "1999"))
alpha_mcd <- as.numeric(Sys.getenv("ALPHA_MCD", "0.75"))
set.seed(as.integer(Sys.getenv("SEED_REAL", "20260927")))

source("R/gv_methods.R")

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required.")
}
if (!requireNamespace("rrcov", quietly = TRUE)) {
  stop("Package 'rrcov' is required.")
}

dat <- MASS::crabs
dat$group <- paste0(dat$sp, dat$sex)
groups <- c("BM", "BF", "OM")
vars <- c("FL", "RW", "CL", "CW", "BD")

if (!all(groups %in% dat$group)) stop("Required crab groups are missing.")
X_list <- lapply(groups, function(g) {
  as.matrix(dat[dat$group == g, vars, drop = FALSE])
})
names(X_list) <- groups

## Fixed, reproducible rowwise contamination:
## 10% of the BM group, shifted in the same direction by 10 group-specific
## standard deviations in every measurement.
X_cont <- X_list
n_cont <- floor(0.10 * nrow(X_cont[["BM"]]))
cont_idx <- sample.int(nrow(X_cont[["BM"]]), n_cont)
shift <- 10 * apply(X_cont[["BM"]], 2, sd)
X_cont[["BM"]][cont_idx, ] <- sweep(
  X_cont[["BM"]][cont_idx, , drop = FALSE], 2, shift, "+"
)

run_one <- function(X, data_label) {
  ans <- run_three_tests(X, B = B, alpha = alpha_mcd)
  data.frame(
    data = data_label,
    method = c("LRT", "BMLRT", "MCD-permutation"),
    statistic = c(ans$lrt$statistic,
                  ans$bmlrt$statistic,
                  ans$mcd_perm$statistic),
    p_value = c(ans$lrt$p.value,
                ans$bmlrt$p.value,
                ans$mcd_perm$p.value),
    decision_alpha_005 = c(ans$lrt$p.value < 0.05,
                           ans$bmlrt$p.value < 0.05,
                           ans$mcd_perm$p.value < 0.05),
    mcd_successful = c(NA_integer_, NA_integer_, ans$mcd_perm$successful),
    mcd_failed = c(NA_integer_, NA_integer_, ans$mcd_perm$failed),
    stringsAsFactors = FALSE
  )
}

results <- rbind(
  run_one(X_list, "Original BM-BF-OM"),
  run_one(X_cont, "10% rowwise-contaminated BM-BF-OM")
)

write.csv(results, out_file, row.names = FALSE)

cat("Wrote:", out_file, "\n")
cat("Groups:", paste(groups, collapse = ", "),
    " | variables:", paste(vars, collapse = ", "), "\n")
cat("B:", B, " | MCD alpha:", alpha_mcd,
    " | contaminated BM rows:", paste(cont_idx, collapse = ","), "\n")
print(results, row.names = FALSE)

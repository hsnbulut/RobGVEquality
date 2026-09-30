args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args)) args[1] else "results/pilot"
output_csv <- if (length(args) > 1) args[2] else "results/summary.csv"

files <- list.files(input_dir, pattern = "^task_.*\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task RDS files found")
raw <- do.call(rbind, lapply(files, readRDS))
raw$alpha <- 0.05
raw$reject_lrt <- as.integer(raw$lrt_p < raw$alpha)
raw$reject_bmlrt <- as.integer(raw$bmlrt_p < raw$alpha)
raw$reject_mcd_perm <- as.integer(raw$mcd_perm_p < raw$alpha)

# aggregate() drops rows whose grouping variable is NA.  Use an explicit
# label so clean scenarios remain in the summary table.
if ("contam_group" %in% names(raw)) {
  raw$contam_group <- ifelse(is.na(raw$contam_group), 0L,
                             as.integer(raw$contam_group))
}

summary <- aggregate(
  cbind(reject_lrt, reject_bmlrt, reject_mcd_perm) ~
    study + k + p + n_i + distribution + shape_difference + delta +
    contamination + epsilon + magnitude + contam_group,
  data = raw, FUN = function(z) mean(z, na.rm = TRUE))
dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, output_csv, row.names = FALSE)
saveRDS(list(raw = raw, summary = summary),
        sub("\\.csv$", ".rds", output_csv))
message("Wrote ", output_csv)

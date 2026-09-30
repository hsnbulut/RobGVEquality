source("R/gv_methods.R")
source("R/simulate_data.R")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}

manifest_path <- get_arg("manifest", "manifest_pilot.rds")
task <- as.integer(get_arg("task", "1"))
reps <- as.integer(get_arg("reps", "20"))
B <- as.integer(get_arg("B", "99"))
alpha <- as.numeric(get_arg("alpha", "0.75"))
out_dir <- get_arg("out", "results/pilot")
base_seed <- as.integer(get_arg("seed", "20260923"))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
manifest <- readRDS(manifest_path)
if (task < 1 || task > nrow(manifest)) stop("task outside manifest")
scenario <- as.list(manifest[task, , drop = FALSE])

res <- vector("list", reps)
for (r in seq_len(reps)) {
  set.seed(base_seed + 100000L * task + r)
  X_list <- simulate_scenario(scenario)
  ans <- tryCatch(run_three_tests(X_list, B = B, alpha = alpha),
                  error = function(e) list(error = conditionMessage(e)))
  res[[r]] <- data.frame(
    task = task, rep = r, study = scenario$study,
    k = scenario$k, p = scenario$p, n_i = scenario$n_i,
    distribution = scenario$distribution,
    shape_difference = scenario$shape_difference,
    delta = scenario$delta, contamination = scenario$contamination,
    epsilon = scenario$epsilon, magnitude = scenario$magnitude,
    contam_group = if (!is.null(scenario$contam_group))
      scenario$contam_group else NA_integer_,
    lrt_p = if (!is.null(ans$lrt)) ans$lrt$p.value else NA_real_,
    bmlrt_p = if (!is.null(ans$bmlrt)) ans$bmlrt$p.value else NA_real_,
    mcd_perm_p = if (!is.null(ans$mcd_perm)) ans$mcd_perm$p.value else NA_real_,
    lrt_stat = if (!is.null(ans$lrt)) ans$lrt$statistic else NA_real_,
    bmlrt_stat = if (!is.null(ans$bmlrt)) ans$bmlrt$statistic else NA_real_,
    mcd_perm_stat = if (!is.null(ans$mcd_perm)) ans$mcd_perm$statistic else NA_real_,
    mcd_perm_successful = if (!is.null(ans$mcd_perm))
      ans$mcd_perm$successful else NA_integer_,
    mcd_perm_failed = if (!is.null(ans$mcd_perm))
      ans$mcd_perm$failed else NA_integer_,
    error = if (!is.null(ans$error)) ans$error else NA_character_,
    stringsAsFactors = FALSE
  )
}

outfile <- file.path(out_dir, sprintf("task_%04d.rds", task))
saveRDS(do.call(rbind, res), outfile)
message("Wrote ", outfile)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "manifest_sensitivity.rds"

# All configurations remain low/moderate dimensional and satisfy p < n_i.
# The sensitivity grid varies sample size, dimension, contamination fraction,
# and signal strength while retaining rowwise contamination only.
make_grid <- function(study, configs, n_i, delta, epsilon,
                      shape_difference, contamination, magnitude = 10) {
  g <- expand.grid(config_id = seq_len(nrow(configs)),
                   n_i = n_i, delta = delta, epsilon = epsilon,
                   KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  data.frame(
    study = study,
    k = configs$k[g$config_id],
    p = configs$p[g$config_id],
    n_i = g$n_i,
    distribution = "normal",
    shape_difference = shape_difference,
    delta = g$delta,
    contamination = contamination,
    epsilon = g$epsilon,
    magnitude = magnitude,
    effect_group = 2L,
    contam_group = if (contamination == "row") 1L else NA_integer_,
    stringsAsFactors = FALSE
  )
}

configs_k2 <- data.frame(k = c(2, 2), p = c(2, 5))
configs_all <- data.frame(k = c(2, 2, 3), p = c(2, 5, 5))

type1_clean <- make_grid(
  "type1_clean", configs_k2, c(30, 60, 100), 1, 0,
  TRUE, "clean", magnitude = 6
)

type1_row <- make_grid(
  "type1_contaminated", configs_all, c(30, 60, 100), 1,
  c(.05, .10, .20), TRUE, "row", magnitude = 10
)

power_clean <- make_grid(
  "power_clean", configs_k2, c(30, 60, 100), 2, 0,
  FALSE, "clean", magnitude = 6
)

# Signal-strength sensitivity at the smallest baseline configuration.
power_delta <- make_grid(
  "power_clean_delta", data.frame(k = 2, p = 2), 30,
  c(1.5, 3), 0, FALSE, "clean", magnitude = 6
)

# Rowwise-contaminated power under the main alternative delta = 2.
power_row <- make_grid(
  "power_contaminated", configs_all, c(30, 60, 100), 2,
  .10, FALSE, "row", magnitude = 10
)

manifest <- rbind(type1_clean, type1_row, power_clean,
                  power_delta, power_row)
manifest$scenario_id <- seq_len(nrow(manifest))
manifest <- manifest[, c("scenario_id", setdiff(names(manifest), "scenario_id"))]

saveRDS(manifest, out)
message("Wrote ", nrow(manifest), " scenarios to ", out)

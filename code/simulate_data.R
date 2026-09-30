shape_matrix <- function(p, id = 1L) {
  if (id == 1L) return(diag(p))
  S <- diag(p)
  r <- if (id == 2L) 4 else 9
  S[1, 1] <- r
  S[2, 2] <- 1 / r
  S
}

make_sigma <- function(p, group_id, gv_root = 1, shape = TRUE) {
  gv_root * shape_matrix(p, if (shape) group_id else 1L)
}

generate_group <- function(n, Sigma, distribution = "normal",
                           contamination = "clean", epsilon = 0,
                           magnitude = 6, contaminate = TRUE) {
  p <- ncol(Sigma)
  Z <- matrix(rnorm(n * p), nrow = n)
  if (distribution %in% c("t3", "t5")) {
    df <- as.numeric(sub("t", "", distribution))
    Z <- Z / sqrt(rchisq(n, df) / df)
  }
  X <- Z %*% chol(Sigma)
  m <- floor(epsilon * n)
  if (m > 0 && contamination != "clean" && contaminate) {
    s <- sqrt(diag(Sigma))
    if (contamination == "row") {
      ii <- sample.int(n, m)
      # One-sided, same-direction contamination creates a high-leverage
      # group-specific distortion of the classical covariance determinant.
      shifts <- rep(1, m)
      X[ii, ] <- X[ii, , drop = FALSE] +
        magnitude * shifts %o% s
    } else if (contamination == "cell") {
      # Use distinct cells and the same sign; sampling with replacement and
      # random signs made the old contamination too weak and non-directional.
      ij <- sample.int(n * p, m * p, replace = FALSE)
      rr <- ((ij - 1L) %% n) + 1L
      cc <- ((ij - 1L) %/% n) + 1L
      jj <- cbind(rr, cc)
      X[jj] <- X[jj] + magnitude *
        s[cc]
    }
  }
  X
}

simulate_scenario <- function(scenario) {
  k <- scenario$k
  n_i <- rep(scenario$n_i, k)
  contam_group <- if (!is.null(scenario$contam_group) &&
                      !is.na(scenario$contam_group)) {
    as.integer(scenario$contam_group)
  } else {
    seq_len(k)
  }
  lapply(seq_len(k), function(i) {
    root <- if (i == scenario$effect_group) scenario$delta else 1
    Sigma <- make_sigma(scenario$p, i, root,
                        shape = scenario$shape_difference)
    generate_group(n_i[i], Sigma, scenario$distribution,
                   scenario$contamination, scenario$epsilon,
                   scenario$magnitude,
                   contaminate = i %in% contam_group)
  })
}

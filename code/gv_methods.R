logdet_spd <- function(S) {
  S <- (S + t(S)) / 2
  R <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(R)) return(NA_real_)
  2 * sum(log(diag(R)))
}

classical_scatter <- function(X) {
  Xc <- sweep(X, 2, colMeans(X), "-")
  crossprod(Xc) / nrow(X)
}

gv_statistic <- function(scatter_list, n_i) {
  p <- nrow(scatter_list[[1]])
  log_a <- vapply(scatter_list, function(S) logdet_spd(S) / p, numeric(1))
  if (any(!is.finite(log_a))) return(NA_real_)
  a <- exp(log_a)
  n <- sum(n_i)
  p * (n * log(sum(n_i * a) / n) - sum(n_i * log(a)))
}

lrt_test <- function(X_list) {
  n_i <- vapply(X_list, nrow, integer(1))
  S <- lapply(X_list, classical_scatter)
  T <- gv_statistic(S, n_i)
  list(statistic = T,
       p.value = pchisq(T, df = length(X_list) - 1L, lower.tail = FALSE))
}

lgamma_multivariate <- function(a, p) {
  p * (p - 1) / 4 * log(pi) +
    sum(lgamma(a + (1 - seq_len(p)) / 2))
}

wishart_det_moment <- function(nu, p, exponent) {
  exponent * p * log(2) +
    lgamma_multivariate(nu / 2 + exponent, p) -
    lgamma_multivariate(nu / 2, p)
}

composition_matrix <- function(total, parts) {
  if (parts == 1L) return(matrix(total, nrow = 1L))
  do.call(rbind, lapply(0:total, function(first) {
    cbind(first, composition_matrix(total - first, parts - 1L))
  }))
}

nth_z_moment <- function(r, j, n_i, p) {
  k <- length(n_i)
  nu <- n_i - 1
  numerator <- setdiff(seq_len(k), j)
  comps <- composition_matrix(r, length(numerator))
  terms <- apply(comps, 1, function(comp) {
    log_coef <- lgamma(r + 1) - sum(lgamma(comp + 1))
    log_prod <- sum(vapply(seq_along(numerator), function(h) {
      wishart_det_moment(nu[numerator[h]], p, comp[h] / p)
    }, numeric(1)))
    exp(log_coef + log_prod)
  })
  denominator <- exp(wishart_det_moment(nu[j], p, -r / p))
  denominator * sum(terms)
}

najarzadeh_phi <- function(n_i, p) {
  k <- length(n_i)
  n <- sum(n_i)
  elog <- numeric(k)
  for (j in seq_len(k)) {
    m <- vapply(1:4, nth_z_moment, numeric(1), j = j, n_i = n_i, p = p)
    mu <- m[1]
    variance <- m[2] - mu^2
    cm3 <- m[3] - 3 * m[2] * mu + 2 * mu^3
    cm4 <- m[4] - 4 * m[3] * mu + 6 * m[2] * mu^2 - 3 * mu^4
    elog[j] <- log1p(mu) - variance / (2 * (1 + mu)^2) +
      2 * cm3 / (6 * (1 + mu)^3) -
      6 * cm4 / (24 * (1 + mu)^4)
  }
  minus_two_log_C <- -n * p * log(n) + p * sum(n_i * log(n_i))
  expected_T <- minus_two_log_C + p * sum(n_i * elog)
  expected_T / (k - 1)
}

bmlrt_test <- function(X_list) {
  n_i <- vapply(X_list, nrow, integer(1))
  p <- ncol(X_list[[1]])
  raw <- lrt_test(X_list)$statistic
  phi <- najarzadeh_phi(n_i, p)
  corrected <- phi * raw
  list(statistic = corrected, raw.statistic = raw, phi = phi,
       p.value = pchisq(corrected, df = length(X_list) - 1L,
                        lower.tail = FALSE))
}

mcd_fit <- function(X, alpha = 0.75) {
  fit <- rrcov::CovMcd(X, alpha = alpha)
  S <- (fit@cov + t(fit@cov)) / 2
  ld <- logdet_spd(S)
  if (!is.finite(ld)) stop("MCD scatter is not positive definite")
  a <- exp(ld / ncol(X))
  list(center = as.numeric(fit@center), scatter = S,
       a = a, shape = S / a)
}

matrix_inv_sqrt <- function(S) {
  ee <- eigen((S + t(S)) / 2, symmetric = TRUE)
  if (any(ee$values <= 0)) stop("shape matrix is not positive definite")
  ee$vectors %*% diag(1 / sqrt(ee$values)) %*% t(ee$vectors)
}

robust_permutation_test <- function(X_list, B = 999, alpha = 0.75) {
  n_i <- vapply(X_list, nrow, integer(1))
  observed_fits <- lapply(X_list, mcd_fit, alpha = alpha)

  # Remove group-specific shape/orientation while preserving the robust
  # generalized-variance scale a_i. Do not divide by sqrt(a_i): that would
  # normalize away the scale differences under the alternative.
  shape_adjusted <- lapply(seq_along(X_list), function(i) {
    fit <- observed_fits[[i]]
    Z <- sweep(X_list[[i]], 2, fit$center, "-")
    Z %*% matrix_inv_sqrt(fit$shape)
  })

  # Refit the observed transformed groups so that T_obs and T_perm use the
  # same statistic and the same transformation pipeline.
  observed_adjusted_fits <- lapply(shape_adjusted, mcd_fit, alpha = alpha)
  T_obs <- gv_statistic(
    lapply(observed_adjusted_fits, function(z) z$scatter), n_i
  )

  pooled <- do.call(rbind, shape_adjusted)
  labels <- rep(seq_along(n_i), n_i)
  T_perm <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    permuted <- pooled[sample.int(nrow(pooled)), , drop = FALSE]
    perm_groups <- lapply(seq_along(n_i), function(i) {
      permuted[which(labels == i), , drop = FALSE]
    })
    fits <- tryCatch(lapply(perm_groups, mcd_fit, alpha = alpha),
                     error = function(e) NULL)
    if (!is.null(fits)) {
      T_perm[b] <- gv_statistic(lapply(fits, function(z) z$scatter), n_i)
    }
  }
  ok <- is.finite(T_perm)
  if (!any(ok)) stop("all MCD permutation fits failed")
  list(statistic = T_obs,
       p.value = (1 + sum(T_perm[ok] >= T_obs)) / (1 + sum(ok)),
       successful = sum(ok), failed = sum(!ok),
       permutation.statistics = T_perm)
}

run_three_tests <- function(X_list, B = 999, alpha = 0.75) {
  list(lrt = lrt_test(X_list),
       bmlrt = bmlrt_test(X_list),
       mcd_perm = robust_permutation_test(X_list, B = B, alpha = alpha))
}

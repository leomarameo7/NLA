# Nested-Library Analysis (Huang, Chang & Hsieh 2024) with GP-EDM instead
# of the original S-map.
#
# - nla_left() / nla_right(): the paper's Algorithm 1 / Algorithm 2.
# - Two-stage fit (Munch & Rogers): hyperparameters are fit once per outer
#   library `n`, then reused (fixedpars=) for every shrinking-library
#   refit -- skips the slow hyperparameter search, keeps only a fast
#   matrix inversion.
# - Lags are built ONCE over the full series (build_lag_frame()): a
#   library point can legitimately use history from before its own start.
# - Scaling is fixed ONCE over the full series, not recomputed per subset
#   (see build_lag_frame()) -- otherwise fixedpars silently stop matching
#   the data they're applied to as the library shrinks.

suppressPackageStartupMessages({
  library(GPEDM)
})

#' Gaussian-kernel smoothing of an error curve
#'
#' Nadaraya-Watson smoother, standing in for the "Gaussian filter" of Huang
#' et al. 2024 Algorithm 1/2 Step 2. `bandwidth` is in units of the index
#' spacing (i.e. of `Dskipstep`).
smooth_gaussian <- function(idx, y, bandwidth) {
  keep <- is.finite(y)
  if (sum(keep) < 3) return(rep(NA_real_, length(y)))
  out <- rep(NA_real_, length(y))
  xi <- idx[keep]; yi <- y[keep]
  for (i in seq_along(idx)) {
    w <- exp(-0.5 * ((xi - idx[i]) / bandwidth)^2)
    out[i] <- sum(w * yi) / sum(w)
  }
  out
}

#' Decide whether a (smoothed) error curve is valley-shaped and, if so,
#' return the location of the minimum.
#'
#' Our discriminant (the original paper's exact "Sec 2.2" discriminant is not
#' reproduced in its supplement): collapse the sign sequence of the first
#' difference of the smoothed curve (ignoring near-zero, numerically flat,
#' differences), and require it to be a single run of "-" followed by a
#' single run of "+" (i.e. exactly one interior local minimum, monotone
#' decrease then monotone increase). The minimum may not sit at either
#' boundary, matching "if it exists" in the algorithm's stated result.
find_valley <- function(idx, y_smooth, flat_tol = 1e-8) {
  ok <- is.finite(y_smooth)
  if (sum(ok) < 3) return(list(is_valley = FALSE, tau_hat = NA_real_))
  idx <- idx[ok]; y_smooth <- y_smooth[ok]
  o <- order(idx); idx <- idx[o]; y_smooth <- y_smooth[o]

  d <- diff(y_smooth)
  scale <- max(abs(y_smooth), flat_tol)
  s <- ifelse(abs(d) < flat_tol * scale, 0, sign(d))
  s_nz <- s[s != 0]
  if (length(s_nz) < 2) return(list(is_valley = FALSE, tau_hat = NA_real_))

  runs <- rle(s_nz)$values
  is_valley <- length(runs) == 2 && identical(runs, c(-1, 1))

  imin <- which.min(y_smooth)
  at_boundary <- imin == 1 || imin == length(y_smooth)
  is_valley <- is_valley && !at_boundary

  list(is_valley = is_valley, tau_hat = if (is_valley) idx[imin] else NA_real_)
}

#' Build lag/predictor columns once over the full series, and z-score them
#' ONCE using a fixed reference (mean/sd over all complete rows).
#'
#' Hyperparameters fixed via `fixedpars` are only meaningful if the scale of
#' the data doesn't change between the fit that produced them and later
#' refits that reuse them. `fitGP(..., scaling="global")` recomputes
#' mean/sd from whatever subset of rows it is given, so as the library
#' shrinks across the inner `l` loop the effective scaling silently drifts
#' out from under the fixed hyperparameters. We instead scale once, up
#' front, over the full series, and use `scaling="none"` in every fitGP()
#' call below so a given library's fixedpars stay valid regardless of how
#' much of the library is later dropped.
#'
#' @return list(lagdf, xnames) where lagdf has a `time` column (1..T) plus
#'   the (scaled) response column `y` and (scaled) lag columns `xnames`.
build_lag_frame <- function(y, E, tau) {
  T <- length(y)
  base <- data.frame(time = seq_len(T), y = y)
  lagdf <- makelags(data = base, y = "y", E = E, tau = tau, append = TRUE)
  xnames <- setdiff(names(lagdf), c("time", "y"))

  complete <- stats::complete.cases(lagdf[, c("y", xnames)])
  ymean <- mean(lagdf$y[complete]); ysd <- sd(lagdf$y[complete])
  lagdf$y <- (lagdf$y - ymean) / ysd
  for (xn in xnames) {
    xmean <- mean(lagdf[[xn]][complete]); xsd <- sd(lagdf[[xn]][complete])
    lagdf[[xn]] <- (lagdf[[xn]] - xmean) / xsd
  }
  list(lagdf = lagdf, xnames = xnames)
}

#' Fit GP hyperparameters on library rows [lib_start, lib_end] (inclusive,
#' in `time` units) and evaluate out-of-sample RMSE on test rows.
fit_hyperparams <- function(lagdf, xnames, lib_start, lib_end, test_rows) {
  train <- lagdf[lagdf$time >= lib_start & lagdf$time <= lib_end, ]
  fit <- tryCatch(
    fitGP(data = train, y = "y", x = xnames, scaling = "none", newdata = test_rows),
    error = function(e) NULL
  )
  fit
}

#' Refit with hyperparameters fixed (Munch/Rogers 2-stage trick): only the
#' matrix inversion for the new (smaller) library is recomputed.
refit_fixed <- function(lagdf, xnames, lib_start, lib_end, test_rows, fixedpars) {
  train <- lagdf[lagdf$time >= lib_start & lagdf$time <= lib_end, ]
  if (nrow(na.omit(train[, c("y", xnames)])) < length(xnames) + 2) return(NA_real_)
  fit <- tryCatch(
    fitGP(data = train, y = "y", x = xnames, scaling = "none",
          fixedpars = fixedpars, newdata = test_rows),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  unname(fit$outsampfitstats["rmse"])
}

#' Nested-Library Analysis, checking the LEFT part of the time series
#' (Algorithm 1 of Huang et al. 2024), with GP-EDM in place of S-map.
#'
#' @param y        Full univariate time series (numeric vector, time-ordered).
#' @param L        Index (1-based) marking the end of the library / start of
#'                  the test set: library = y[1:L], test = y[(L+1):R].
#' @param R        Last index used (defaults to length(y)).
#' @param E, tau    Embedding dimension and lag for GP-EDM.
#' @param Dskipstep Step size for both the outer (n) and inner (l) loops.
#' @param bandwidth Gaussian-smoothing bandwidth (in index units) applied to
#'                  each error curve before valley detection.
#' @return list(tau_hat = median change-point estimate (NA if none found),
#'              tau_n = per-outer-n estimates, curves = list of per-n error
#'              curves for diagnostic plotting)
nla_left <- function(y, L, R = length(y), E = 3, tau = 1, Dskipstep = 5,
                      bandwidth = 2 * Dskipstep, min_lib = E * tau + 3) {
  lf <- build_lag_frame(y, E, tau)
  lagdf <- lf$lagdf; xnames <- lf$xnames
  test_rows <- lagdf[lagdf$time >= (L + 1) & lagdf$time <= R, ]

  n_seq <- seq(0, L - min_lib, by = Dskipstep)
  tau_n <- rep(NA_real_, length(n_seq))
  curves <- vector("list", length(n_seq))

  for (i in seq_along(n_seq)) {
    n <- n_seq[i]
    lib_start_n <- max(1, n)
    hp <- fit_hyperparams(lagdf, xnames, lib_start_n, L, test_rows)
    if (is.null(hp)) next
    fixedpars <- hp$pars[seq_len(length(xnames) + 2)]

    l_seq <- seq(n, L - min_lib, by = Dskipstep)
    if (length(l_seq) < 3) next
    err <- vapply(l_seq, function(l) {
      refit_fixed(lagdf, xnames, max(1, l), L, test_rows, fixedpars)
    }, numeric(1))

    y_smooth <- smooth_gaussian(l_seq, err, bandwidth)
    v <- find_valley(l_seq, y_smooth)
    tau_n[i] <- v$tau_hat
    curves[[i]] <- data.frame(n = n, l = l_seq, rmse = err, rmse_smooth = y_smooth)
  }

  list(tau_hat = if (all(is.na(tau_n))) NA_real_ else median(tau_n, na.rm = TRUE),
       tau_n = data.frame(n = n_seq, tau_hat = tau_n),
       curves = do.call(rbind, curves))
}

#' Nested-Library Analysis, checking the RIGHT part of the time series
#' (Algorithm 2 of Huang et al. 2024), with GP-EDM in place of S-map.
#'
#' @param L Index marking the end of the test set / start of the library:
#'          test = y[1:L], library = y[(L+1):R].
nla_right <- function(y, L, R = length(y), E = 3, tau = 1, Dskipstep = 5,
                       bandwidth = 2 * Dskipstep, min_lib = E * tau + 3) {
  lf <- build_lag_frame(y, E, tau)
  lagdf <- lf$lagdf; xnames <- lf$xnames
  test_rows <- lagdf[lagdf$time >= 1 & lagdf$time <= L, ]

  n_seq <- seq(R, L + 1 + min_lib, by = -Dskipstep)
  tau_n <- rep(NA_real_, length(n_seq))
  curves <- vector("list", length(n_seq))

  for (i in seq_along(n_seq)) {
    n <- n_seq[i]
    hp <- fit_hyperparams(lagdf, xnames, L + 1, n, test_rows)
    if (is.null(hp)) next
    fixedpars <- hp$pars[seq_len(length(xnames) + 2)]

    l_seq <- seq(n, L + 1 + min_lib, by = -Dskipstep)
    if (length(l_seq) < 3) next
    err <- vapply(l_seq, function(l) {
      refit_fixed(lagdf, xnames, L + 1, l, test_rows, fixedpars)
    }, numeric(1))

    y_smooth <- smooth_gaussian(l_seq, err, bandwidth)
    v <- find_valley(l_seq, y_smooth)
    tau_n[i] <- v$tau_hat
    curves[[i]] <- data.frame(n = n, l = l_seq, rmse = err, rmse_smooth = y_smooth)
  }

  list(tau_hat = if (all(is.na(tau_n))) NA_real_ else median(tau_n, na.rm = TRUE),
       tau_n = data.frame(n = n_seq, tau_hat = tau_n),
       curves = do.call(rbind, curves))
}

#' Run both the left- and right-checking NLA passes and combine them.
#'
#' Mirrors how Huang et al. 2024 use two symmetric test sets (one at each
#' end of the series) as a robustness check (their S8 Text): `L_left` sets
#' library/test split for the left-checking algorithm (library = 1:L_left,
#' test = (L_left+1):R_end), and `L_right` sets the split for the
#' right-checking algorithm (test = 1:L_right, library = (L_right+1):R_end).
nla_gpedm <- function(y, L_left, L_right, R_end = length(y),
                       E = 3, tau = 1, Dskipstep = 5, bandwidth = 2 * Dskipstep) {
  left <- nla_left(y, L = L_left, R = R_end, E = E, tau = tau,
                    Dskipstep = Dskipstep, bandwidth = bandwidth)
  right <- nla_right(y, L = L_right, R = R_end, E = E, tau = tau,
                       Dskipstep = Dskipstep, bandwidth = bandwidth)
  ests <- c(left$tau_hat, right$tau_hat)
  list(tau_hat = if (all(is.na(ests))) NA_real_ else median(ests, na.rm = TRUE),
       left = left, right = right)
}

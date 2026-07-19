# Synthetic regime-shift test system for NLA-GPEDM.
#
# Huang, Chang & Hsieh (2024, PLOS Comp Biol, "Detecting shifts in nonlinear
# dynamics using EDM with Nested-Library Analysis") validate the original
# (S-map based) NLA algorithm on a 4-variable forced food-chain model (their
# Table A, S3 Text): a resource N, forced by a slowly increasing driver
# a(t) through a Holling type-III ("hill function") feedback that makes the
# N subsystem bistable over an intermediate range of a, feeding a 3-species
# chain x -> y -> z whose vital rates depend on N. Ramping a(t) across the
# upper fold of the N subsystem produces one abrupt regime shift that
# propagates through the whole food chain.
#
# We tried to reproduce that exact 4-variable model from the OCR'd
# supplementary PDF, but at least one of its fine-grained rate constants is
# not reliably recoverable from the scan: as OCR'd, the forcing ramp
# c4=0.375 sweeps a from 0.25 to >1800 over the paper's own t in [0,5000]
# (instantly blowing past any bifurcation), and with d1 = c0 N + c1
# (c0=0.75, c1=3.7) the loss rate d1 (~3.9) always exceeds f1's own upper
# bound (~0.13), so y collapses to 0 and x diverges for every initial
# condition (verified numerically -- see git history of this file). We
# can't tell which digits were corrupted without the original source, so we
# do not present invented replacement constants as "the paper's parameters".
#
# We keep the part of the model we CAN take directly, unambiguously, from
# the table -- the forced bistable N subsystem, which is the actual source
# of the regime shift in the full model:
#
#   dN/dt = a(t) - b N + c N^m/(N^m+1)      b=0.8, m=8, c=1  (as OCR'd)
#   a(t)  = c4 t + c5                        c5=0.25 as OCR'd
#
# (c4, the ramp slope, is chosen by us -- see find_fold_points() -- so that
# a(t) sweeps across the subsystem's upper fold once, at a chosen point in
# the sampled window: we derive the one constant that wasn't recoverable
# from the OCR'd text, instead of guessing digits.)
#
# For the downstream food chain we use the well-established, numerically
# stable Hastings & Powell (1991) 3-species model instead of Huang et al.'s
# (numerically inconsistent, as extracted) x/y/z equations. This is the
# same model bundled as `HastPow3sp` in the GPEDM package itself (same
# a=5,b=3,c=0.1,d=2,m=0.4,mu=0.01 parameterisation, documented there as
# chaotic), and is the standard EDM-community benchmark for "a control
# parameter drives a food chain between fixed-point / cyclic / chaotic
# regimes" -- precisely the kind of test Medeiros et al. 2025 (PNAS) used
# GP-EDM on. We drive the attack-rate parameter `a` of the Hastings-Powell
# chain with the (rescaled) trajectory of the bistable N subsystem, so that
# N's abrupt fold-jump forces an abrupt shift in the food chain's dynamical
# regime, at a precisely known time:
#
#   dX/dt = X(1-X) - a(t) XY/(1+b X)
#   dY/dt = a(t) XY/(1+b X) - c YZ/(1+d Y) - m Y
#   dZ/dt = c YZ/(1+d Y) - mu Z
#   a(t)  = a_lo + (a_hi - a_lo) * clip((N_raw(t)-N_lo)/(N_hi-N_lo), 0, 1)
#
# where N_raw(t) is the (noisy) bistable-subsystem trajectory and N_lo/N_hi
# are its pre-/post-shift branch levels. Poisson-arrival multiplicative
# process shocks and additive Gaussian measurement noise are added to all
# state variables, as described in Huang et al. 2024 S3 Text.

suppressPackageStartupMessages({
  library(deSolve)
})

bistable_rhs <- function(t, state, params) {
  N <- state["N"]
  with(as.list(params), {
    a <- c4 * t + c5
    dN <- a - b * N + c * N^m / (N^m + 1)
    list(c(N = dN))
  })
}

#' Locate the saddle-node (fold) points of dN/dt = a - b N + c N^m/(N^m+1),
#' i.e. the values of forcing `a` at which the number of equilibria changes.
find_fold_points <- function(b = 0.8, m = 8, c = 1, Nmax = 3, n = 20000) {
  N <- seq(1e-3, Nmax, length.out = n)
  g <- b * N - c * N^m / (N^m + 1)   # equilibrium condition: a = g(N)
  d <- diff(g)
  turns <- which(diff(sign(d)) != 0) + 1
  data.frame(N = N[turns], a = g[turns])
}

default_bistable_params <- function(c4 = 0.000743, c5 = 0.25) {
  list(b = 0.8, m = 8, c = 1, c4 = c4, c5 = c5)
}

default_hp_params <- function(a_lo = 2.2, a_hi = 5.2) {
  list(a_lo = a_lo, a_hi = a_hi, b = 3, c = 0.1, d = 2, m = 0.4, mu = 0.01)
}

#' Simulate the forced bistable N subsystem alone (process noise only, no
#' measurement noise -- used internally to build the a(t) forcing driving
#' the food chain, and to define the ground-truth change point).
simulate_bistable <- function(Tmax, dt_sample, lambda, rho,
                               state0 = c(N = 0.30),
                               params = default_bistable_params(),
                               seed) {
  set.seed(seed)
  n_events <- rpois(1, lambda * Tmax)
  event_times <- sort(runif(n_events, 1e-6, Tmax - 1e-6))
  shocks <- matrix(runif(n_events, 1 - rho, 1 + rho), ncol = 1)

  event_func <- function(t, state, params) {
    i <- which(abs(event_times - t) < 1e-6)[1]
    if (!is.na(i)) state <- pmax(state * shocks[i, ], 1e-6)
    state
  }

  samp_times <- seq(0, Tmax, by = dt_sample)
  sol <- ode(y = state0, times = samp_times, func = bistable_rhs, parms = params,
             method = "ode45",
             events = if (n_events > 0) list(func = event_func, time = event_times) else NULL)
  traj <- as.data.frame(sol)
  traj[complete.cases(traj), ]
}

hp_rhs_forced <- function(t, state, params) {
  with(as.list(c(state, params)), {
    a <- a_fun(t)
    dX <- X * (1 - X) - a * X * Y / (1 + b * X)
    dY <- a * X * Y / (1 + b * X) - c * Y * Z / (1 + d * Y) - m * Y
    dZ <- c * Y * Z / (1 + d * Y) - mu * Z
    list(c(X = dX, Y = dY, Z = dZ))
  })
}

#' Simulate the full synthetic regime-shift system: a bistable N subsystem
#' (forced by a slowly ramping control parameter, producing one abrupt
#' fold-jump) driving the attack-rate parameter of a Hastings-Powell
#' 3-species food chain, which therefore undergoes one abrupt shift in
#' dynamical regime at the same, precisely known time.
#'
#' @param Tmax        Total integration time.
#' @param dt_sample   Observation interval (state saved every dt_sample).
#' @param lambda      Poisson rate of process-error shocks (per time unit).
#' @param rho         Half-width of the multiplicative process shock.
#' @param meas_sd      SD of additive Gaussian measurement noise (applied
#'                     after independently z-scoring each variable).
#' @param bistable_params  Params for the N subsystem (see
#'                     default_bistable_params).
#' @param hp_params    Params for the food chain (see default_hp_params).
#' @param state0_hp    Initial state for the food chain, c(X=,Y=,Z=).
#' @param branch_frac  Fraction of the series (at each end) used to
#'                     estimate the N subsystem's pre-/post-shift branch
#'                     levels, for rescaling N into the food chain's a(t).
#' @param seed         RNG seed for reproducibility.
simulate_foodchain <- function(Tmax = 900, dt_sample = 2, lambda = 0.01,
                                rho = 0.08, meas_sd = 0.03,
                                bistable_params = default_bistable_params(),
                                hp_params = default_hp_params(),
                                state0_hp = c(X = 0.5, Y = 0.1, Z = 9),
                                branch_frac = 0.15,
                                seed = 1) {
  n_traj <- simulate_bistable(Tmax, dt_sample, lambda, rho,
                               params = bistable_params, seed = seed)

  n_pts <- nrow(n_traj)
  lo_idx <- seq_len(floor(branch_frac * n_pts))
  hi_idx <- seq(ceiling((1 - branch_frac) * n_pts), n_pts)
  N_lo <- mean(n_traj$N[lo_idx])
  N_hi <- mean(n_traj$N[hi_idx])

  N_interp <- approxfun(n_traj$time, n_traj$N, rule = 2)
  a_fun <- function(t) {
    frac <- (N_interp(t) - N_lo) / (N_hi - N_lo)
    hp_params$a_lo + (hp_params$a_hi - hp_params$a_lo) * pmin(pmax(frac, 0), 1)
  }

  set.seed(seed + 10000)
  n_events <- rpois(1, lambda * Tmax)
  event_times <- sort(runif(n_events, 1e-6, Tmax - 1e-6))
  shocks <- matrix(runif(3 * n_events, 1 - rho, 1 + rho), ncol = 3)
  event_func <- function(t, state, params) {
    i <- which(abs(event_times - t) < 1e-6)[1]
    if (!is.na(i)) state <- pmax(state * shocks[i, ], 1e-6)
    state
  }

  samp_times <- seq(0, Tmax, by = dt_sample)
  params <- c(hp_params, list(a_fun = a_fun))
  sol <- ode(y = state0_hp, times = samp_times, func = hp_rhs_forced, parms = params,
             method = "ode45",
             events = if (n_events > 0) list(func = event_func, time = event_times) else NULL)
  traj <- as.data.frame(sol)
  traj <- traj[complete.cases(traj), ]

  traj$N_raw <- N_interp(traj$time)
  traj$a <- a_fun(traj$time)

  for (v in c("X", "Y", "Z")) {
    traj[[v]] <- scale(traj[[v]])[, 1] + rnorm(nrow(traj), 0, meas_sd)
  }
  rownames(traj) <- NULL
  traj
}

#' Empirically determine the regime-shift time index from a simulated
#' trajectory, as "the moment when the absolute rate of change of N reaches
#' its maximum" (Huang et al. 2024 S3 Text), returned as a 1-based row index.
detect_true_changepoint <- function(traj, var = "N_raw") {
  d <- abs(diff(traj[[var]]))
  which.max(d) + 1L
}

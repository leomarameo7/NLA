# Synthetic regime-shift benchmark for NLA-GPEDM.
#
# What: a 2-species predator-prey model (Rosenzweig-MacArthur) whose
#       carrying capacity K is slowly ramped up, crossing the classic
#       "paradox of enrichment" Hopf bifurcation -> one clean shift from a
#       stable equilibrium to a limit cycle, at a precisely known time.
# Why:  simplest possible system with a textbook, unambiguous regime
#       shift -- good for checking NLA-GPEDM recovers a KNOWN answer.
#
#   dV/dt = r V (1 - V/K) - a V P / (1 + a h V)
#   dP/dt = e a V P / (1 + a h V) - m P
#   K(t)  = K_lo + (K_hi - K_lo) * t / Tmax
#
# With r=1, a=1.1, h=0.8, e=0.7, m=0.4 (standard textbook values), the
# coexistence equilibrium is stable for K < K_crit = 3.0 and unstable
# (-> limit cycle) for K > K_crit (found by bisection on simulated
# oscillation amplitude). Ramping K from 2.0 to 4.0 crosses K_crit exactly
# halfway through the simulation.

suppressPackageStartupMessages({
  library(deSolve)
})

rm_rhs <- function(t, state, params) {
  with(as.list(c(state, params)), {
    K <- K_lo + (K_hi - K_lo) * t / Tmax
    dV <- r * V * (1 - V / K) - a * V * P / (1 + a * h * V)
    dP <- e * a * V * P / (1 + a * h * V) - m * P
    list(c(V = dV, P = dP))
  })
}

default_rm_params <- function(K_lo = 2, K_hi = 4, Tmax = 500) {
  list(r = 1, a = 1.1, h = 0.8, e = 0.7, m = 0.4,
       K_lo = K_lo, K_hi = K_hi, Tmax = Tmax)
}

#' Simulate the forced Rosenzweig-MacArthur model, with Poisson-arrival
#' multiplicative process shocks and additive Gaussian measurement noise.
#'
#' @param Tmax, dt_sample  Total time and sampling interval.
#' @param lambda, rho      Poisson shock rate and multiplicative half-width.
#' @param meas_sd          SD of measurement noise (after z-scoring).
#' @param params           See default_rm_params().
#' @param seed             RNG seed.
simulate_regime_shift <- function(Tmax = 500, dt_sample = 1, lambda = 0.01,
                                   rho = 0.05, meas_sd = 0.02,
                                   params = default_rm_params(Tmax = Tmax),
                                   seed = 1) {
  set.seed(seed)
  n_events <- rpois(1, lambda * Tmax)
  event_times <- sort(runif(n_events, 1e-6, Tmax - 1e-6))
  shocks <- matrix(runif(2 * n_events, 1 - rho, 1 + rho), ncol = 2)
  event_func <- function(t, state, params) {
    i <- which(abs(event_times - t) < 1e-6)[1]
    if (!is.na(i)) state <- pmax(state * shocks[i, ], 1e-6)
    state
  }

  samp_times <- seq(0, Tmax, by = dt_sample)
  sol <- ode(y = c(V = 1, P = 0.5), times = samp_times, func = rm_rhs, parms = params,
             method = "ode45",
             events = if (n_events > 0) list(func = event_func, time = event_times) else NULL)
  traj <- as.data.frame(sol)
  traj <- traj[complete.cases(traj), ]

  traj$K <- with(params, K_lo + (K_hi - K_lo) * traj$time / Tmax)
  for (v in c("V", "P")) {
    traj[[v]] <- scale(traj[[v]])[, 1] + rnorm(nrow(traj), 0, meas_sd)
  }
  rownames(traj) <- NULL
  traj
}

#' True change-point index: the (1-based) row where K(t) crosses K_crit.
true_changepoint <- function(traj, K_crit = 3.0) {
  which.min(abs(traj$K - K_crit))
}

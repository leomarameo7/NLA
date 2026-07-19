# Validate NLA-GPEDM on the synthetic regime-shift system (see
# simulate_foodchain.R): a bistable resource N, forced by a slowly ramping
# control parameter, drives the attack-rate parameter of a Hastings-Powell
# 3-species food chain from a stable equilibrium into a limit cycle at a
# precisely known time. This mirrors the validation Huang et al. 2024 did
# for the original S-map based NLA, and the "control parameter drives a
# regime shift" framing Medeiros et al. 2025 used GP-EDM for.
#
# Run from the repository root: Rscript R_code/run_synthetic_test.R

suppressPackageStartupMessages({
  library(here)
})
source(here("R_code", "simulate_foodchain.R"))
source(here("R_code", "nla_gpedm.R"))

dir.create(here("results"), showWarnings = FALSE)

traj <- simulate_foodchain(seed = 1)
cp_true <- detect_true_changepoint(traj)
T <- nrow(traj)
cat(sprintf("Simulated series length: %d, true change point (index): %d\n", T, cp_true))

varname <- "Y"
y <- traj[[varname]]

E <- 3; tau <- 1; Dskipstep <- 8; bandwidth <- 16
L_left <- 300
L_right <- 80

cat("Running left-checking NLA-GPEDM (Algorithm 1)...\n")
t0 <- Sys.time()
left <- nla_left(y, L = L_left, R = T, E = E, tau = tau,
                  Dskipstep = Dskipstep, bandwidth = bandwidth)
cat(sprintf("  done in %.1fs, tau_hat = %s\n", as.numeric(Sys.time() - t0, units = "secs"), left$tau_hat))

cat("Running right-checking NLA-GPEDM (Algorithm 2)...\n")
t0 <- Sys.time()
right <- nla_right(y, L = L_right, R = T, E = E, tau = tau,
                     Dskipstep = Dskipstep, bandwidth = bandwidth)
cat(sprintf("  done in %.1fs, tau_hat = %s\n", as.numeric(Sys.time() - t0, units = "secs"), right$tau_hat))

combined <- median(c(left$tau_hat, right$tau_hat), na.rm = TRUE)

cat(sprintf("\nTrue change point:      %d\n", cp_true))
cat(sprintf("Left-check estimate:    %s\n", left$tau_hat))
cat(sprintf("Right-check estimate:   %s\n", right$tau_hat))
cat(sprintf("Combined estimate:      %s\n", combined))

saveRDS(list(traj = traj, cp_true = cp_true, varname = varname,
             left = left, right = right, combined = combined,
             params = list(E = E, tau = tau, Dskipstep = Dskipstep,
                            bandwidth = bandwidth, L_left = L_left, L_right = L_right)),
        here("results", "synthetic_nla_gpedm.rds"))

cat("\nSaved results to results/synthetic_nla_gpedm.rds\n")

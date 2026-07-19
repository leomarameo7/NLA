# Validate NLA-GPEDM against a KNOWN change point: a 2-species
# predator-prey model whose carrying capacity is ramped across a Hopf
# bifurcation (paradox of enrichment), shifting from a stable equilibrium
# to a limit cycle (see simulate_regime_shift.R).
#
# Run from the repository root: Rscript R_code/run_synthetic_test.R

suppressPackageStartupMessages({
  library(here)
})
source(here("R_code", "simulate_regime_shift.R"))
source(here("R_code", "nla_gpedm.R"))

dir.create(here("results"), showWarnings = FALSE)

traj <- simulate_regime_shift(seed = 1)
cp_true <- true_changepoint(traj)
T <- nrow(traj)
cat(sprintf("Series length: %d, true change point (index): %d\n", T, cp_true))

y <- traj$V
E <- 3; tau <- 1; Dskipstep <- 8; bandwidth <- 16
L_left <- 350
L_right <- 150

cat("Left-checking NLA-GPEDM...\n")
t0 <- Sys.time()
left <- nla_left(y, L = L_left, R = T, E = E, tau = tau, Dskipstep = Dskipstep, bandwidth = bandwidth)
cat(sprintf("  done in %.1fs, tau_hat = %s\n", as.numeric(Sys.time() - t0, units = "secs"), left$tau_hat))

cat("Right-checking NLA-GPEDM...\n")
t0 <- Sys.time()
right <- nla_right(y, L = L_right, R = T, E = E, tau = tau, Dskipstep = Dskipstep, bandwidth = bandwidth)
cat(sprintf("  done in %.1fs, tau_hat = %s\n", as.numeric(Sys.time() - t0, units = "secs"), right$tau_hat))

combined <- median(c(left$tau_hat, right$tau_hat), na.rm = TRUE)

cat(sprintf("\nTrue change point:    %d\n", cp_true))
cat(sprintf("Left-check estimate:  %s\n", left$tau_hat))
cat(sprintf("Right-check estimate: %s\n", right$tau_hat))
cat(sprintf("Combined estimate:    %s\n", combined))

saveRDS(list(traj = traj, cp_true = cp_true, left = left, right = right, combined = combined,
             params = list(E = E, tau = tau, Dskipstep = Dskipstep, bandwidth = bandwidth,
                            L_left = L_left, L_right = L_right)),
        here("results", "synthetic_nla_gpedm.rds"))

cat("\nSaved results to results/synthetic_nla_gpedm.rds\n")

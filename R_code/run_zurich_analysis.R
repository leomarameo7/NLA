# Apply NLA-GPEDM to the Lake Zurich monthly plankton time series (Merz et
# al. 2023) and compare the detected regime-shift time against the
# change-point analysis of Medeiros et al. 2025 (PNAS), who found:
#   - Large green algae (Gr2): November 1987 (their most well-resolved shift)
#   - Small cryptophytes (Cr1): October 1988
#   - Omnivores (Om):           June 1986
# suggesting an ecosystem-wide shift around 1986-1988.
#
# Run from the repository root: Rscript R_code/run_zurich_analysis.R

suppressPackageStartupMessages({
  library(here)
})
source(here("R_code", "process_zurich_data.R"))
source(here("R_code", "nla_gpedm.R"))

dir.create(here("results"), showWarnings = FALSE)

z <- process_zurich_data(here("data", "raw", "lake_zurich.csv"),
                          start_year = 1977, end_year = 1999)
T <- nrow(z)
cat(sprintf("Processed series: %d months (%s to %s)\n", T,
            format(min(z$date)), format(max(z$date))))

focal_vars <- list(
  Gr2 = list(col = "log_Gr2", ref_year = 1987, ref_month = 11, label = "Green algae (large)"),
  Cr1 = list(col = "log_Cr1", ref_year = 1988, ref_month = 10, label = "Cryptophytes (small)"),
  Om  = list(col = "log_Om",  ref_year = 1986, ref_month = 6,  label = "Omnivores")
)

E <- 3; tau <- 1; Dskipstep <- 8; bandwidth <- 16
L_left <- 200
L_right <- 60

results_list <- list()
for (nm in names(focal_vars)) {
  spec <- focal_vars[[nm]]
  y <- z[[spec$col]]
  ref_idx <- zurich_month_index(z, spec$ref_year, spec$ref_month)

  cat(sprintf("\n=== %s (%s) ===\n", nm, spec$label))
  cat(sprintf("Medeiros et al. 2025 reference change point: %d-%02d (index %d)\n",
              spec$ref_year, spec$ref_month, ref_idx))

  t0 <- Sys.time()
  left <- nla_left(y, L = L_left, R = T, E = E, tau = tau,
                    Dskipstep = Dskipstep, bandwidth = bandwidth)
  right <- nla_right(y, L = L_right, R = T, E = E, tau = tau,
                       Dskipstep = Dskipstep, bandwidth = bandwidth)
  combined <- median(c(left$tau_hat, right$tau_hat), na.rm = TRUE)
  dt <- as.numeric(Sys.time() - t0, units = "secs")

  cat(sprintf("Left-check: %s | Right-check: %s | Combined: %s | (%.1fs)\n",
              left$tau_hat, right$tau_hat, combined, dt))

  results_list[[nm]] <- list(varname = nm, label = spec$label, col = spec$col,
                              ref_idx = ref_idx, ref_date = sprintf("%d-%02d", spec$ref_year, spec$ref_month),
                              left = left, right = right, combined = combined)
}

saveRDS(list(z = z, results = results_list,
             params = list(E = E, tau = tau, Dskipstep = Dskipstep,
                            bandwidth = bandwidth, L_left = L_left, L_right = L_right)),
        here("results", "zurich_nla_gpedm.rds"))

cat("\nSaved results to results/zurich_nla_gpedm.rds\n")

summary_df <- do.call(rbind, lapply(results_list, function(r) {
  data.frame(variable = r$label, medeiros_date = r$ref_date, medeiros_index = r$ref_idx,
             nla_gpedm_left = r$left$tau_hat, nla_gpedm_right = r$right$tau_hat,
             nla_gpedm_combined = r$combined)
}))
print(summary_df, row.names = FALSE)
write.csv(summary_df, here("results", "zurich_nla_gpedm_summary.csv"), row.names = FALSE)

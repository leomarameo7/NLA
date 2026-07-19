# Process the Lake Zurich monthly plankton time series (Merz et al. 2023)
# for NLA-GPEDM analysis, following the variable naming and (log-transform)
# convention used by Medeiros et al. 2025 for the same dataset.

suppressPackageStartupMessages({
  library(zoo)
})

zurich_long_names <- c(
  Cy = "Cyanobacteria", Gr1 = "Green algae (small)", Gr2 = "Green algae (large)",
  Di1 = "Diatoms (small)", Di2 = "Diatoms (large)", Go1 = "Gold algae (small)",
  Go2 = "Gold algae (large)", Cr1 = "Cryptophytes (small)", Cr2 = "Cryptophytes (large)",
  Mi = "Mixotrophic flagellates", He = "Large herbivores", Om = "Omnivores",
  Pr = "Invertebrate predators", temperature = "Temperature", phosphorus = "Phosphate"
)

#' Load and process the raw Lake Zurich CSV: builds a monthly time index,
#' linearly interpolates the handful of scattered single-month gaps, and
#' (optionally) restricts to a date window.
#'
#' @param path        Path to lake_zurich.csv.
#' @param start_year  First year to keep (inclusive), or NULL for all data.
#' @param end_year    Last year to keep (inclusive), or NULL for all data.
process_zurich_data <- function(path, start_year = NULL, end_year = NULL) {
  z <- read.csv(path)
  z <- z[order(z$year, z$month), ]
  z$time_idx <- seq_len(nrow(z))
  z$date <- as.Date(sprintf("%d-%02d-01", z$year, z$month))

  plankton_vars <- setdiff(names(zurich_long_names), c("temperature", "phosphorus"))
  for (v in c(plankton_vars, "temperature", "phosphorus")) {
    z[[v]] <- na.approx(z[[v]], x = z$time_idx, na.rm = FALSE)
  }
  z <- z[complete.cases(z[, plankton_vars]), ]

  for (v in plankton_vars) {
    z[[paste0("log_", v)]] <- log1p(z[[v]])
  }

  if (!is.null(start_year)) z <- z[z$year >= start_year, ]
  if (!is.null(end_year)) z <- z[z$year <= end_year, ]
  z$time_idx <- seq_len(nrow(z))
  rownames(z) <- NULL
  z
}

#' Row index (1-based, within the processed/filtered data frame `z`) closest
#' to a given year-month, for marking reference change-point dates on plots.
zurich_month_index <- function(z, year, month) {
  target <- as.Date(sprintf("%d-%02d-01", year, month))
  which.min(abs(z$date - target))
}

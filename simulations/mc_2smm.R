#!/usr/bin/env Rscript
# =============================================================================
# Monte Carlo Simulation for 2SMM Validation
# =============================================================================
#
# Validates the 2SMM method (Wall & Amemiya 2000, 2003) across:
#   - Sample sizes: 200, 500, 1000, 2000
#   - Interaction degrees: 2 (X:Z), 3 (X:X:Z), 4 (X:X:X:Z), complex (multi-term)
#   - Error distribution assumption: "general" vs "normal"
#   - Data normality: normal vs non-normal latent/error distributions
#
# Total: 4 x 4 x 2 x 2 = 64 conditions x 500 reps = 32,000 model fits
#
# Usage:
#   Rscript simulations/mc_2smm.R
#
# Output:
#   simulations/results/cond_XX.rds   (per-condition raw results)
#   simulations/results/mc_2smm_summary.csv  (aggregated metrics)
# =============================================================================

cat("=============================================================\n")
cat("  Monte Carlo Simulation: 2SMM Validation\n")
cat("=============================================================\n\n")

# ---- Setup ----
# Use devtools::load_all() to load the development version of modsem
# (modsem_2smm is not yet on CRAN). Falls back to library(modsem) if
# devtools is unavailable or the source tree is not found.
if (requireNamespace("devtools", quietly = TRUE) &&
    file.exists(file.path("DESCRIPTION"))) {
  devtools::load_all(".", quiet = TRUE)
} else {
  library(modsem)
}

set.seed(42)

RESULTS_DIR <- file.path("simulations", "results")
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

N_REPS <- 500L
CI_Z <- qnorm(0.975)  # 1.96 for 95% CI

# Number of cores for parallelization
N_CORES <- if (.Platform$OS.type == "windows") {
  1L
} else {
  4L
}
cat(sprintf("Using %d cores for parallel execution\n\n", N_CORES))


# =============================================================================
# DGP Parameters
# =============================================================================

# Common measurement model: 3 indicators per factor, all loadings = 1
# Measurement error variance = 0.16 (sd = 0.4) per indicator
LOADING <- 1.0
INTERCEPT <- 0.0
ERROR_VAR <- 0.16

# Latent covariance: Cov(X, Z) = 0.2, Var(X) = Var(Z) = 1
COV_XZ <- 0.2

# Structural error variance
ZETA_VAR <- 0.25

# True parameter values per model type
TRUE_PARAMS <- list(
  degree2 = list(
    syntax = "Y ~ X + Z + X:Z",
    beta = c("Y~X" = 0.6, "Y~Z" = 0.4, "Y~X:Z" = 0.3),
    dgp = function(X, Z) 0.6 * X + 0.4 * Z + 0.3 * X * Z
  ),
  degree3 = list(
    syntax = "Y ~ X + Z + X:X:Z",
    beta = c("Y~X" = 0.5, "Y~Z" = 0.3, "Y~X:X:Z" = 0.3),
    dgp = function(X, Z) 0.5 * X + 0.3 * Z + 0.3 * X^2 * Z
  ),
  degree4 = list(
    syntax = "Y ~ X + Z + X:X:X:Z",
    beta = c("Y~X" = 0.5, "Y~Z" = 0.3, "Y~X:X:X:Z" = 0.3),
    dgp = function(X, Z) 0.5 * X + 0.3 * Z + 0.3 * X^3 * Z
  ),
  complex = list(
    syntax = "Y ~ X + Z + X:Z + X:X + Z:Z + X:Z:Z",
    beta = c("Y~X" = 0.5, "Y~Z" = 0.3, "Y~X:Z" = 0.3,
             "Y~X:X" = 0.2, "Y~Z:Z" = 0.15, "Y~X:Z:Z" = 0.2),
    dgp = function(X, Z) {
      0.5 * X + 0.3 * Z + 0.3 * X * Z + 0.2 * X^2 + 0.15 * Z^2 + 0.2 * X * Z^2
    }
  )
)


# =============================================================================
# Condition grid
# =============================================================================

conditions <- expand.grid(
  N = c(200L, 500L, 1000L, 2000L),
  degree = c("degree2", "degree3", "degree4", "complex"),
  error.dist = c("general", "normal"),
  data.normal = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
conditions$cond_id <- seq_len(nrow(conditions))

cat(sprintf("Total conditions: %d\n", nrow(conditions)))
cat(sprintf("Reps per condition: %d\n", N_REPS))
cat(sprintf("Total model fits: %d\n\n", nrow(conditions) * N_REPS))


# =============================================================================
# Data generation function
# =============================================================================

generate_data <- function(N, degree, data.normal, seed) {
  set.seed(seed)

  model_info <- TRUE_PARAMS[[degree]]

  # ---- Generate latent variables X, Z ----
  if (data.normal) {
    # Bivariate normal: Var(X)=1, Var(Z)=1, Cov(X,Z)=0.2
    Sigma_LV <- matrix(c(1, COV_XZ, COV_XZ, 1), nrow = 2)
    LV <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_LV)
    X <- LV[, 1]
    Z <- LV[, 2]
  } else {
    # Non-normal: chi-squared(5), centered and scaled to mean=0, var=1
    # Then induce correlation via Cholesky
    X_raw <- (rchisq(N, df = 5) - 5) / sqrt(2 * 5)
    Z_raw <- (rchisq(N, df = 5) - 5) / sqrt(2 * 5)
    # Cholesky to induce Cov(X,Z) = 0.2
    L <- chol(matrix(c(1, COV_XZ, COV_XZ, 1), nrow = 2))
    LV <- cbind(X_raw, Z_raw) %*% L
    X <- LV[, 1]
    Z <- LV[, 2]
  }

  # ---- Structural equation: Y = f(X, Z) + zeta ----
  Y_star <- model_info$dgp(X, Z)
  zeta <- rnorm(N, mean = 0, sd = sqrt(ZETA_VAR))
  Y <- Y_star + zeta

  # ---- Measurement model: 3 indicators per factor ----
  gen_indicators <- function(LV_scores, prefix, data.normal) {
    n <- length(LV_scores)
    if (data.normal) {
      errors <- matrix(rnorm(n * 3, mean = 0, sd = sqrt(ERROR_VAR)), nrow = n, ncol = 3)
    } else {
      # t(5) scaled to variance = 0.16
      # Var(t_5) = 5/(5-2) = 5/3, so scale by sqrt(0.16 / (5/3)) = sqrt(0.096)
      errors <- matrix(rt(n * 3, df = 5) * sqrt(ERROR_VAR / (5 / 3)),
                        nrow = n, ncol = 3)
    }
    indicators <- LOADING * LV_scores + INTERCEPT + errors
    colnames(indicators) <- paste0(prefix, 1:3)
    indicators
  }

  x_ind <- gen_indicators(X, "x", data.normal)
  z_ind <- gen_indicators(Z, "z", data.normal)
  y_ind <- gen_indicators(Y, "y", data.normal)

  as.data.frame(cbind(x_ind, z_ind, y_ind))
}


# =============================================================================
# Single replication function
# =============================================================================

run_one_rep <- function(cond_row, rep_id) {
  N <- cond_row$N
  degree <- cond_row$degree
  error.dist <- cond_row$error.dist
  data.normal <- cond_row$data.normal
  cond_id <- cond_row$cond_id

  # Deterministic seed: condition_id * 10000 + rep_id
  seed <- cond_id * 10000L + rep_id

  model_info <- TRUE_PARAMS[[degree]]
  true_beta <- model_info$beta
  param_names <- names(true_beta)

  # Full model syntax
  model_syntax <- paste0(
    "\n",
    "  X =~ x1 + x2 + x3\n",
    "  Z =~ z1 + z2 + z3\n",
    "  Y =~ y1 + y2 + y3\n",
    "  ", model_info$syntax, "\n"
  )

  result <- tryCatch({
    dat <- generate_data(N, degree, data.normal, seed)

    fit <- modsem_2smm(model_syntax, data = dat,
                       error.dist = error.dist, verbose = FALSE)

    # Extract estimates and SEs from parameter table
    pt <- fit$parTable
    struct <- pt[pt$op == "~", ]

    # Match parameter names: "Y~X", "Y~X:Z", etc.
    est_names <- paste0(struct$lhs, "~", struct$rhs)
    est_vals <- struct$est
    se_vals <- struct$std.error
    names(est_vals) <- est_names
    names(se_vals) <- est_names

    # Build result row for each parameter
    rows <- lapply(param_names, function(pname) {
      est <- if (pname %in% names(est_vals)) est_vals[[pname]] else NA_real_
      se  <- if (pname %in% names(se_vals))  se_vals[[pname]]  else NA_real_
      data.frame(
        cond_id = cond_id, rep_id = rep_id, N = N, degree = degree,
        error.dist = error.dist, data.normal = data.normal,
        param = pname, true_value = true_beta[[pname]],
        est = est, se = se, converged = TRUE,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, rows)

  }, error = function(e) {
    # Return NA rows on failure
    rows <- lapply(param_names, function(pname) {
      data.frame(
        cond_id = cond_id, rep_id = rep_id, N = N, degree = degree,
        error.dist = error.dist, data.normal = data.normal,
        param = pname, true_value = true_beta[[pname]],
        est = NA_real_, se = NA_real_, converged = FALSE,
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })

  result
}


# =============================================================================
# Main simulation loop
# =============================================================================

cat("Starting simulation...\n")
cat("=============================================================\n\n")

total_start <- proc.time()

for (i in seq_len(nrow(conditions))) {
  cond <- conditions[i, ]
  result_file <- file.path(RESULTS_DIR, sprintf("cond_%02d.rds", i))

  # Skip if already completed (resume support)
  if (file.exists(result_file)) {
    cat(sprintf("[%2d/%d] SKIP (already exists): N=%4d, %s, error.dist=%s, normal=%s\n",
                i, nrow(conditions), cond$N, cond$degree, cond$error.dist, cond$data.normal))
    next
  }

  cat(sprintf("[%2d/%d] Running: N=%4d, %s, error.dist=%s, normal=%s ... ",
              i, nrow(conditions), cond$N, cond$degree, cond$error.dist, cond$data.normal))

  cond_start <- proc.time()

  # Run replications in parallel (or sequentially on Windows)
  rep_results <- if (N_CORES > 1L) {
    parallel::mclapply(seq_len(N_REPS), function(rep_id) {
      run_one_rep(cond, rep_id)
    }, mc.cores = N_CORES)
  } else {
    lapply(seq_len(N_REPS), function(rep_id) {
      run_one_rep(cond, rep_id)
    })
  }

  # Combine into single data frame
  cond_results <- do.call(rbind, rep_results)

  # Save incrementally
  saveRDS(cond_results, result_file)

  elapsed <- (proc.time() - cond_start)[["elapsed"]]
  n_converged <- sum(cond_results$converged[!duplicated(cond_results$rep_id)])
  cat(sprintf("done (%.1fs, %d/%d converged)\n", elapsed, n_converged, N_REPS))
}

total_elapsed <- (proc.time() - total_start)[["elapsed"]]
cat(sprintf("\nAll conditions complete. Total time: %.1f minutes\n\n",
            total_elapsed / 60))


# =============================================================================
# Summary: compute metrics
# =============================================================================

cat("=============================================================\n")
cat("  Computing summary metrics\n")
cat("=============================================================\n\n")

# Read all condition results
all_files <- sort(list.files(RESULTS_DIR, pattern = "^cond_\\d+\\.rds$", full.names = TRUE))
all_results <- do.call(rbind, lapply(all_files, readRDS))

# Compute metrics per condition x parameter
summary_list <- list()

cond_param_groups <- split(all_results,
                           list(all_results$cond_id, all_results$param),
                           drop = TRUE)

for (grp in cond_param_groups) {
  cond_id <- grp$cond_id[1]
  param <- grp$param[1]
  true_val <- grp$true_value[1]

  # Filter to converged reps
  conv <- grp[grp$converged, ]
  n_total <- nrow(grp)
  n_converged <- nrow(conv)
  convergence_rate <- n_converged / n_total

  if (n_converged < 2) {
    summary_list[[length(summary_list) + 1]] <- data.frame(
      cond_id = cond_id, N = grp$N[1], degree = grp$degree[1],
      error.dist = grp$error.dist[1], data.normal = grp$data.normal[1],
      param = param, true_value = true_val,
      bias = NA, rel_bias_pct = NA, empirical_sd = NA,
      mean_se = NA, se_sd_ratio = NA, rmse = NA, coverage_95 = NA,
      convergence_rate = convergence_rate, n_converged = n_converged,
      stringsAsFactors = FALSE
    )
    next
  }

  est <- conv$est
  se <- conv$se

  # Remove any NA estimates (extra safety)
  valid <- !is.na(est) & !is.na(se)
  est <- est[valid]
  se <- se[valid]
  n_valid <- length(est)

  if (n_valid < 2) next

  bias <- mean(est) - true_val
  rel_bias_pct <- (bias / true_val) * 100
  empirical_sd <- sd(est)
  mean_se <- mean(se)
  se_sd_ratio <- mean_se / empirical_sd
  rmse <- sqrt(mean((est - true_val)^2))

  # 95% CI coverage
  ci_lower <- est - CI_Z * se
  ci_upper <- est + CI_Z * se
  coverage_95 <- mean(true_val >= ci_lower & true_val <= ci_upper)

  summary_list[[length(summary_list) + 1]] <- data.frame(
    cond_id = cond_id, N = grp$N[1], degree = grp$degree[1],
    error.dist = grp$error.dist[1], data.normal = grp$data.normal[1],
    param = param, true_value = true_val,
    bias = bias, rel_bias_pct = rel_bias_pct, empirical_sd = empirical_sd,
    mean_se = mean_se, se_sd_ratio = se_sd_ratio, rmse = rmse,
    coverage_95 = coverage_95,
    convergence_rate = convergence_rate, n_converged = n_converged,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)
rownames(summary_df) <- NULL

# Save summary
summary_file <- file.path(RESULTS_DIR, "mc_2smm_summary.csv")
write.csv(summary_df, summary_file, row.names = FALSE)
cat(sprintf("Summary saved to: %s\n\n", summary_file))


# =============================================================================
# Print formatted results
# =============================================================================

cat("=============================================================\n")
cat("  Results Summary\n")
cat("=============================================================\n\n")

# Print by degree
for (deg in c("degree2", "degree3", "degree4", "complex")) {
  sub <- summary_df[summary_df$degree == deg, ]
  if (nrow(sub) == 0) next

  cat(sprintf("--- %s ---\n", toupper(deg)))
  cat(sprintf("%-10s %-8s %-8s %-6s %7s %8s %7s %7s %8s %7s %8s %5s\n",
              "Param", "N", "err.dist", "normal",
              "Bias", "RelB(%)", "EmpSD", "MeanSE", "SE/SD", "RMSE", "Cov95", "Conv%"))
  cat(paste(rep("-", 105), collapse = ""), "\n")

  for (j in seq_len(nrow(sub))) {
    r <- sub[j, ]
    cat(sprintf("%-10s %-8d %-8s %-6s %7.4f %8.2f %7.4f %7.4f %8.3f %7.4f %8.3f %5.1f\n",
                r$param, r$N, r$error.dist,
                ifelse(r$data.normal, "Y", "N"),
                r$bias, r$rel_bias_pct, r$empirical_sd, r$mean_se,
                r$se_sd_ratio, r$rmse, r$coverage_95,
                r$convergence_rate * 100))
  }
  cat("\n")
}


# =============================================================================
# Print key diagnostics
# =============================================================================

cat("=============================================================\n")
cat("  Key Diagnostics\n")
cat("=============================================================\n\n")

# Convergence rates by degree
cat("Convergence rates by degree:\n")
conv_by_deg <- tapply(summary_df$convergence_rate, summary_df$degree, mean, na.rm = TRUE)
for (d in names(conv_by_deg)) {
  cat(sprintf("  %s: %.1f%%\n", d, conv_by_deg[[d]] * 100))
}

# SE/SD ratio summary
cat("\nSE/SD ratio (should be ~1.0):\n")
sesd_by_N <- tapply(summary_df$se_sd_ratio, summary_df$N, median, na.rm = TRUE)
for (n in names(sesd_by_N)) {
  cat(sprintf("  N=%s: median SE/SD = %.3f\n", n, sesd_by_N[[n]]))
}

# Coverage summary
cat("\n95% CI coverage (should be ~0.95):\n")
cov_by_N <- tapply(summary_df$coverage_95, summary_df$N, median, na.rm = TRUE)
for (n in names(cov_by_N)) {
  cat(sprintf("  N=%s: median coverage = %.3f\n", n, cov_by_N[[n]]))
}

# Non-normal comparison
cat("\nNon-normal data: general vs normal error.dist (median |bias|):\n")
nn <- summary_df[!summary_df$data.normal, ]
if (nrow(nn) > 0) {
  bias_gen <- median(abs(nn$bias[nn$error.dist == "general"]), na.rm = TRUE)
  bias_norm <- median(abs(nn$bias[nn$error.dist == "normal"]), na.rm = TRUE)
  cat(sprintf("  error.dist='general': median |bias| = %.4f\n", bias_gen))
  cat(sprintf("  error.dist='normal':  median |bias| = %.4f\n", bias_norm))
}

cat("\nDone.\n")

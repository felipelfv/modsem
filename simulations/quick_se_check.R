#!/usr/bin/env Rscript
# Quick validation: does the full sandwich Ω fix SE underestimation?
# Runs degree-2 model (Y ~ X + Z + X:Z) with N=500, 200 reps

cat("Quick SE validation (full sandwich Ω)\n")
cat("======================================\n\n")

devtools::load_all(".", quiet = TRUE)
set.seed(42)

N <- 500L
N_REPS <- 200L
CI_Z <- qnorm(0.975)

LOADING <- 1.0
ERROR_VAR <- 0.16
COV_XZ <- 0.2
ZETA_VAR <- 0.25

TRUE_BETA <- c("Y~X" = 0.6, "Y~Z" = 0.4, "Y~X:Z" = 0.3)

model_syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
"

generate_data <- function(N, seed) {
  set.seed(seed)
  Sigma_LV <- matrix(c(1, COV_XZ, COV_XZ, 1), nrow = 2)
  LV <- MASS::mvrnorm(N, mu = c(0, 0), Sigma = Sigma_LV)
  X <- LV[, 1]; Z <- LV[, 2]
  Y <- 0.6 * X + 0.4 * Z + 0.3 * X * Z + rnorm(N, 0, sqrt(ZETA_VAR))

  gen_ind <- function(lv, prefix) {
    err <- matrix(rnorm(N * 3, 0, sqrt(ERROR_VAR)), N, 3)
    ind <- lv + err
    colnames(ind) <- paste0(prefix, 1:3)
    ind
  }
  as.data.frame(cbind(gen_ind(X, "x"), gen_ind(Z, "z"), gen_ind(Y, "y")))
}

cat(sprintf("N=%d, reps=%d, model=degree2 (Y ~ X + Z + X:Z)\n\n", N, N_REPS))

results <- vector("list", N_REPS)

for (rep_id in seq_len(N_REPS)) {
  if (rep_id %% 20 == 0) cat(sprintf("  rep %d/%d\n", rep_id, N_REPS))

  res <- tryCatch({
    dat <- generate_data(N, seed = 10000L + rep_id)
    fit <- modsem_2smm(model_syntax, data = dat, error.dist = "normal", verbose = FALSE)

    pt <- fit$parTable
    struct <- pt[pt$op == "~", ]
    est_names <- paste0(struct$lhs, "~", struct$rhs)
    est_vals <- setNames(struct$est, est_names)
    se_vals  <- setNames(struct$std.error, est_names)

    data.frame(
      rep_id = rep_id,
      param = names(TRUE_BETA),
      true_value = unname(TRUE_BETA),
      est = unname(est_vals[names(TRUE_BETA)]),
      se  = unname(se_vals[names(TRUE_BETA)]),
      converged = TRUE,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      rep_id = rep_id,
      param = names(TRUE_BETA),
      true_value = unname(TRUE_BETA),
      est = NA_real_, se = NA_real_, converged = FALSE,
      stringsAsFactors = FALSE
    )
  })

  results[[rep_id]] <- res
}

all_res <- do.call(rbind, results)

cat("\n======================================\n")
cat("Results (full sandwich Ω):\n")
cat("======================================\n\n")

cat(sprintf("%-10s %7s %7s %7s %7s %7s %7s\n",
            "Param", "Bias", "EmpSD", "MeanSE", "SE/SD", "RMSE", "Cov95"))
cat(paste(rep("-", 62), collapse = ""), "\n")

for (pname in names(TRUE_BETA)) {
  sub <- all_res[all_res$param == pname & all_res$converged, ]
  if (nrow(sub) < 2) next

  est <- sub$est; se <- sub$se
  valid <- !is.na(est) & !is.na(se)
  est <- est[valid]; se <- se[valid]

  true_val <- TRUE_BETA[[pname]]
  bias <- mean(est) - true_val
  emp_sd <- sd(est)
  mean_se <- mean(se)
  se_sd <- mean_se / emp_sd
  rmse <- sqrt(mean((est - true_val)^2))
  cov95 <- mean(true_val >= est - CI_Z * se & true_val <= est + CI_Z * se)

  cat(sprintf("%-10s %7.4f %7.4f %7.4f %7.3f %7.4f %7.3f\n",
              pname, bias, emp_sd, mean_se, se_sd, rmse, cov95))
}

n_conv <- sum(all_res$converged[!duplicated(all_res$rep_id)])
cat(sprintf("\nConverged: %d/%d (%.1f%%)\n", n_conv, N_REPS, 100 * n_conv / N_REPS))
cat("\nTarget: SE/SD ~ 1.0, Coverage ~ 0.95\n")
cat("(Old Ω* had SE/SD ~ 0.82-0.87, Coverage ~ 0.88-0.91)\n")

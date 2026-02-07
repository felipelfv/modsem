devtools::load_all()

# ===========================================================================
# Mathematical validation of 2SMM bias-corrected SE
# Tests that the code implements Wall & Amemiya (2000, 2003) correctly:
#   - Error moments match Isserlis' theorem under normality
#   - Per-obs bias-corrected terms average to unbiased_moment
#   - Bias-corrected ell has mean zero (algebraic identity)
#   - Full sandwich SE >= Omega*-only SE (PSD property)
#   - Bias-corrected ell differs from raw ell
#   - Monte Carlo calibration: SE/SD ~ 1, coverage ~ 0.95
# ===========================================================================


# ---- Setup: run the internal 2SMM pipeline on oneInt ----

m1 <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
"

pt <- modsemify(m1)
lvs <- getLVs(pt)
etas <- getSortedEtas(pt, isLV = TRUE)
syntaxCFA <- parTableToSyntax(pt[pt$op == "=~", ])

# Stage 1: CFA + Bartlett scores
stage1 <- estimate2SMM_stage1(syntaxCFA, data = oneInt, lvs = lvs)

# Error moments under normality
errMom <- estimateErrorMoments(
  W = stage1$W, Lambda = stage1$Lambda,
  Theta = stage1$Theta, v = stage1$v,
  factor_names = stage1$factor_names,
  max_order = 4, error.dist = "normal"
)
eps_cumulants <- errMom$eps_cumulants
error_moment_fn <- make_error_moment_fn(stage1$W, eps_cumulants)
factor_names <- stage1$factor_names

# T_hat: n * lavaan::vcov (ACOV of CFA params)
param_info <- get_cfa_relevant_params(stage1$cfa_fit)
full_vcov <- lavaan::vcov(stage1$cfa_fit)
T_hat <- stage1$n * full_vcov[param_info$par_idx, param_info$par_idx, drop = FALSE]

# Stage 2: bias-corrected MoM regression
stage2 <- estimate2SMM_stage2(
  f_hat = stage1$f_hat, Sigma_ee = stage1$Sigma_ee,
  eta_name = "Y", linear_preds = c("X", "Z"), int_terms = "X:Z",
  error_moment_fn = error_moment_fn, factor_names = factor_names
)

reg_multiidx <- lapply(
  c("(Intercept)", "X", "Z", "X:Z"),
  term_to_multiindex, factor_names = factor_names
)
names(reg_multiidx) <- c("(Intercept)", "X", "Z", "X:Z")

# Bias-corrected ell (Paper Eq 13: ell_i = m_hat_i - M_hat_i alpha_hat)
bc_ell <- compute_bc_ell(
  stage1$f_hat, stage2$alpha_hat, reg_multiidx,
  "Y", factor_names, error_moment_fn
)

# C_hat via numerical differentiation (Paper Eq 14)
C_hat <- compute_C_hat(
  stage1 = stage1, alpha_hat = stage2$alpha_hat,
  reg_multiidx = reg_multiidx, eta_name = "Y",
  factor_names = factor_names, lvs = lvs,
  param_info = param_info, eps_cumulants = eps_cumulants
)


# ===========================================================================
# Test 1: Error moments under normality (Isserlis' theorem)
#
# For e = W*eps with eps ~ N(0, Theta):
#   Cov(e) = W Theta W' = Sigma_ee = (Lambda' Theta^{-1} Lambda)^{-1}
#   E[e_a] = 0
#   E[e_a^2] = Sigma_ee[a,a]
#   E[e_a e_b] = Sigma_ee[a,b]
#   E[odd-order] = 0
#   E[e_a^4] = 3 Sigma_ee[a,a]^2
#   E[e_a^2 e_b^2] = Sigma_ee[a,a]*Sigma_ee[b,b] + 2*Sigma_ee[a,b]^2
# ===========================================================================

W <- stage1$W
Sigma_ee <- stage1$Sigma_ee
q <- nrow(W)

# E[e_a^2] = Sigma_ee[a,a]
for (a in seq_len(q)) {
  testthat::expect_equal(
    compute_error_joint_moment(c(a, a), W, eps_cumulants),
    Sigma_ee[a, a], tolerance = 1e-10
  )
}

# E[e_a * e_b] = Sigma_ee[a,b] for a != b
for (a in 1:(q - 1)) {
  for (b in (a + 1):q) {
    testthat::expect_equal(
      compute_error_joint_moment(c(a, b), W, eps_cumulants),
      Sigma_ee[a, b], tolerance = 1e-10
    )
  }
}

# E[e_a] = 0
for (a in seq_len(q)) {
  testthat::expect_equal(
    compute_error_joint_moment(a, W, eps_cumulants), 0
  )
}

# E[e_a^3] = 0, E[e_a*e_b*e_c] = 0 (all odd-order moments)
testthat::expect_equal(
  compute_error_joint_moment(c(1L, 1L, 1L), W, eps_cumulants), 0,
  tolerance = 1e-10
)
testthat::expect_equal(
  compute_error_joint_moment(c(1L, 2L, 3L), W, eps_cumulants), 0,
  tolerance = 1e-10
)

# E[e_a^4] = 3 * Sigma_ee[a,a]^2
for (a in seq_len(q)) {
  testthat::expect_equal(
    compute_error_joint_moment(c(a, a, a, a), W, eps_cumulants),
    3 * Sigma_ee[a, a]^2, tolerance = 1e-10
  )
}

# E[e_a^2 * e_b^2] = Sigma_ee[a,a]*Sigma_ee[b,b] + 2*Sigma_ee[a,b]^2
for (a in 1:(q - 1)) {
  for (b in (a + 1):q) {
    testthat::expect_equal(
      compute_error_joint_moment(c(a, a, b, b), W, eps_cumulants),
      Sigma_ee[a, a] * Sigma_ee[b, b] + 2 * Sigma_ee[a, b]^2,
      tolerance = 1e-10
    )
  }
}


# ===========================================================================
# Test 2: mean(compute_bc_individual(r)) == unbiased_moment(r)
#
# This is the key identity: per-observation bias-corrected contributions
# must average to the scalar bias-corrected moment. This ensures
# M_hat = (1/n) sum M_hat_i and m_hat = (1/n) sum m_hat_i.
# ===========================================================================

test_indices <- list(
  c(X = 1L, Z = 0L, Y = 0L),   # order 1 (linear)
  c(X = 0L, Z = 1L, Y = 0L),   # order 1
  c(X = 1L, Z = 1L, Y = 0L),   # order 2 (interaction X*Z)
  c(X = 2L, Z = 0L, Y = 0L),   # order 2 (X^2)
  c(X = 1L, Z = 0L, Y = 1L),   # order 2 (cross X*Y)
  c(X = 0L, Z = 2L, Y = 0L),   # order 2 (Z^2)
  c(X = 2L, Z = 1L, Y = 0L),   # order 3 (X^2*Z)
  c(X = 1L, Z = 1L, Y = 1L),   # order 3 (X*Z*Y)
  c(X = 0L, Z = 2L, Y = 1L)    # order 3 (Z^2*Y)
)

memo_ub <- new.env(hash = TRUE, parent = emptyenv())
for (midx in test_indices) {
  memo_bc <- new.env(hash = TRUE, parent = emptyenv())
  memo_sc <- new.env(hash = TRUE, parent = emptyenv())
  bc_vec <- compute_bc_individual(midx, stage1$f_hat, error_moment_fn, memo_bc, memo_sc)
  ub_val <- unbiased_moment(midx, stage1$f_hat, error_moment_fn, memo_ub)
  testthat::expect_equal(mean(bc_vec), ub_val, tolerance = 1e-12)
}


# ===========================================================================
# Test 3: bias-corrected ell has column means == 0
#
# Since alpha_hat = M_hat^{-1} m_hat, the estimating equation is
# exactly satisfied at the population level:
#   (1/n) sum ell_i = m_hat - M_hat * alpha_hat
#                   = m_hat - M_hat * M_hat^{-1} * m_hat = 0
# ===========================================================================

for (j in seq_len(ncol(bc_ell))) {
  testthat::expect_equal(mean(bc_ell[, j]), 0, tolerance = 1e-10)
}


# ===========================================================================
# Test 4: Full sandwich SE >= Omega*-only SE
#
# Omega = Omega* + C_hat T_hat C_hat'
# Since T_hat is PSD and C_hat T_hat C_hat' is PSD, we have:
#   V_full = (1/n) M^{-1} Omega M^{-1} >= (1/n) M^{-1} Omega* M^{-1}
# elementwise on the diagonal (i.e., all SEs are at least as large)
# ===========================================================================

se_full <- compute2SMMSE(
  alpha_hat = stage2$alpha_hat,
  M_hat_inv = stage2$M_hat_inv,
  R = stage2$R,
  f_eta = stage2$f_eta,
  n = stage1$n,
  C_hat = C_hat,
  T_hat = T_hat,
  ell = bc_ell
)

se_star_only <- compute2SMMSE(
  alpha_hat = stage2$alpha_hat,
  M_hat_inv = stage2$M_hat_inv,
  R = stage2$R,
  f_eta = stage2$f_eta,
  n = stage1$n,
  C_hat = NULL,
  T_hat = NULL,
  ell = bc_ell
)

# Full SE >= Omega*-only SE for every parameter
testthat::expect_true(all(se_full$se >= se_star_only$se - 1e-15))

# C_hat T_hat C_hat' must be PSD (eigenvalues >= 0)
CTCt <- C_hat %*% T_hat %*% t(C_hat)
eig_vals <- eigen(CTCt, symmetric = TRUE, only.values = TRUE)$values
testthat::expect_true(all(eig_vals >= -1e-10))

# The correction should be non-trivial (Stage 1 uncertainty matters)
testthat::expect_true(max(se_full$se / se_star_only$se) > 1.01)


# ===========================================================================
# Test 5: bias-corrected ell differs from raw ell
#
# Raw ell: ell_i = e_hat_i * R_i (residual times regressors)
# BC ell:  ell_i = m_hat_i - M_hat_i alpha_hat
# They must differ because M_hat_i includes bias correction terms
# (subtracting error moments from raw products)
# ===========================================================================

raw_ell <- compute_ell_contributions(
  stage1$f_hat, stage2$alpha_hat, reg_multiidx,
  "Y", factor_names, error_moment_fn = NULL
)

testthat::expect_false(isTRUE(all.equal(bc_ell, raw_ell, tolerance = 1e-6)))


# ===========================================================================
# Test 6: make_error_moment_fn consistency
#
# The factory-created closure should return identical values to
# direct calls to compute_error_joint_moment
# ===========================================================================

emf <- make_error_moment_fn(W, eps_cumulants)

test_s_vecs <- list(
  c(2L, 0L, 0L),
  c(0L, 2L, 0L),
  c(0L, 0L, 2L),
  c(1L, 1L, 0L),
  c(2L, 2L, 0L),
  c(1L, 1L, 1L, 0L, 0L, 0L)  # longer s to test padding
)

for (s in test_s_vecs) {
  # Build idx_vec the same way make_error_moment_fn does internally
  idx_vec <- integer(0)
  for (l in seq_along(s)) {
    if (s[l] > 0) idx_vec <- c(idx_vec, rep(l, s[l]))
  }
  direct_val <- compute_error_joint_moment(idx_vec, W, eps_cumulants)
  fn_val <- emf(s)
  testthat::expect_equal(fn_val, direct_val, tolerance = 1e-15)
}

# Caching: calling twice should return the same result
s_test <- c(1L, 1L, 0L)
val1 <- emf(s_test)
val2 <- emf(s_test)
testthat::expect_identical(val1, val2)


# ===========================================================================
# Test 7: Monte Carlo SE calibration (100 reps, N=400)
#
# Gold-standard validation: if V_hat{alpha} is a correct covariance
# estimator, then:
#   (a) SE/SD ratio ~ 1.0 (mean SE should match empirical SD)
#   (b) 95% CI coverage ~ 0.95
#   (c) Bias should be negligible
#
# Old Omega*-only gave SE/SD ~ 0.82-0.87, Coverage ~ 0.88-0.91.
# The full sandwich should correct this to ~ 1.0 and ~ 0.95.
# ===========================================================================

set.seed(42)
N_SIM   <- 400L
N_REPS  <- 100L
CI_Z    <- qnorm(0.975)

TRUE_BETA <- c("Y~X" = 0.6, "Y~Z" = 0.4, "Y~X:Z" = 0.3)

mc_syntax <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
"

gen_mc_data <- function(N, seed) {
  set.seed(seed)
  Sigma_LV <- matrix(c(1, 0.2, 0.2, 1), 2)
  LV <- MASS::mvrnorm(N, c(0, 0), Sigma_LV)
  X <- LV[, 1]; Z <- LV[, 2]
  Y <- 0.6 * X + 0.4 * Z + 0.3 * X * Z + rnorm(N, 0, sqrt(0.25))
  mk <- function(lv, pf) {
    ind <- lv + matrix(rnorm(N * 3, 0, sqrt(0.16)), N, 3)
    colnames(ind) <- paste0(pf, 1:3)
    ind
  }
  as.data.frame(cbind(mk(X, "x"), mk(Z, "z"), mk(Y, "y")))
}

mc_res <- vector("list", N_REPS)
for (i in seq_len(N_REPS)) {
  mc_res[[i]] <- tryCatch({
    d <- gen_mc_data(N_SIM, 10000L + i)
    fit <- modsem_2smm(mc_syntax, data = d, error.dist = "normal", verbose = FALSE)
    s <- fit$parTable[fit$parTable$op == "~", ]
    nms <- paste0(s$lhs, "~", s$rhs)
    data.frame(
      param = names(TRUE_BETA),
      est = unname(setNames(s$est, nms)[names(TRUE_BETA)]),
      se  = unname(setNames(s$std.error, nms)[names(TRUE_BETA)]),
      ok = TRUE, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(param = names(TRUE_BETA),
               est = NA_real_, se = NA_real_, ok = FALSE,
               stringsAsFactors = FALSE)
  })
}

mc_all <- do.call(rbind, mc_res)

for (pname in names(TRUE_BETA)) {
  sub <- mc_all[mc_all$param == pname & mc_all$ok, ]
  good <- !is.na(sub$est) & !is.na(sub$se)
  est <- sub$est[good]; se <- sub$se[good]
  tv <- TRUE_BETA[[pname]]

  emp_sd  <- sd(est)
  mean_se <- mean(se)
  se_sd   <- mean_se / emp_sd
  cov95   <- mean(tv >= est - CI_Z * se & tv <= est + CI_Z * se)
  bias    <- mean(est) - tv

  # SE/SD ratio should be close to 1.0
  testthat::expect_true(
    se_sd > 0.80 && se_sd < 1.30,
    label = sprintf("%s: SE/SD = %.3f, expected in [0.80, 1.30]", pname, se_sd)
  )

  # 95% CI coverage should be close to 0.95
  testthat::expect_true(
    cov95 >= 0.85,
    label = sprintf("%s: Coverage = %.3f, expected >= 0.85", pname, cov95)
  )

  # Bias should be small relative to sampling variability
  testthat::expect_true(
    abs(bias) < 3 * emp_sd / sqrt(length(est)),
    label = sprintf("%s: |Bias| = %.4f, expected < %.4f",
                    pname, abs(bias), 3 * emp_sd / sqrt(length(est)))
  )
}

# Convergence rate should be high
conv_ok <- vapply(mc_res, function(r) all(r$ok), logical(1))
testthat::expect_true(
  mean(conv_ok) > 0.90,
  label = sprintf("Convergence rate = %.1f%%, expected > 90%%",
                  100 * mean(conv_ok))
)

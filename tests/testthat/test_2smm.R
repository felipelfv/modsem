devtools::load_all()

# ===========================================================================
# Test 1: Basic interaction model via modsem_2smm()
# ===========================================================================
m1 <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z
"

est1 <- modsem_2smm(m1, data = oneInt, verbose = TRUE)
testthat::expect_s3_class(est1, "modsem_2smm")
testthat::expect_s3_class(est1, "modsem")

# ===========================================================================
# Test 2: Dispatch via modsem(..., method = "2smm")
# ===========================================================================
est2 <- modsem(m1, data = oneInt, method = "2smm", verbose = FALSE)
testthat::expect_s3_class(est2, "modsem_2smm")

# Estimates should match direct call
testthat::expect_equal(coef(est1), coef(est2))

# ===========================================================================
# Test 3: S3 methods work
# ===========================================================================

# summary
s <- summary(est1)
testthat::expect_s3_class(s, "summary_modsem_2smm")
print(s)  # should print without error

# coef
co <- coef(est1)
testthat::expect_true(is.numeric(co))
testthat::expect_true(length(co) > 0)
testthat::expect_true("Y~X" %in% names(co))
testthat::expect_true("Y~Z" %in% names(co))
testthat::expect_true("Y~X:Z" %in% names(co))

# vcov
v <- vcov(est1)
testthat::expect_true(is.matrix(v))
testthat::expect_true(nrow(v) > 0)

# nobs
testthat::expect_equal(nobs(est1), nrow(oneInt))

# parameter_estimates
pe <- parameter_estimates(est1)
testthat::expect_true(is.data.frame(pe))
testthat::expect_true("lhs" %in% names(pe))
testthat::expect_true("est" %in% names(pe))

# is_interaction_model
testthat::expect_true(is_interaction_model(est1))

# modsem_inspect
testthat::expect_true(is.data.frame(modsem_inspect(est1, "partable")))
testthat::expect_equal(modsem_inspect(est1, "N"), nrow(oneInt))

# ===========================================================================
# Test 4: Interaction estimate sign agrees with LMS
# ===========================================================================
est_lms <- modsem(m1, data = oneInt, method = "lms", verbose = FALSE)
coef_lms  <- coef(est_lms)
coef_2smm <- coef(est1)

# Interaction effect should have the same sign
int_lms  <- coef_lms[grepl("X:Z", names(coef_lms))][[1]]
int_2smm <- coef_2smm[["Y~X:Z"]]
testthat::expect_equal(sign(int_lms), sign(int_2smm))

# ===========================================================================
# Test 5: error.dist = "general" vs "normal" both work
# ===========================================================================
est_gen  <- modsem_2smm(m1, data = oneInt, error.dist = "general", verbose = FALSE)
est_norm <- modsem_2smm(m1, data = oneInt, error.dist = "normal", verbose = FALSE)

testthat::expect_s3_class(est_gen, "modsem_2smm")
testthat::expect_s3_class(est_norm, "modsem_2smm")

# With normally distributed data, results should be similar
testthat::expect_equal(coef(est_gen)[["Y~X:Z"]],
                       coef(est_norm)[["Y~X:Z"]],
                       tolerance = 0.3)

# ===========================================================================
# Test 6: Backward compatibility snapshot (degree-2)
# ===========================================================================
# Reference values from the original hard-coded implementation
ref_coefs <- c("Y.(Intercept)" = -0.1439088, "Y.X" = 0.6741089,
               "Y.Z" = 0.5610313, "Y.X:Z" = 0.7155527)

# Results must match within numerical tolerance (different code paths)
testthat::expect_equal(unname(est1$coef.structural),
                       unname(ref_coefs),
                       tolerance = 1e-5)

# ===========================================================================
# Test 7: Cubic interaction term (X:X:Z)
# ===========================================================================
m_cubic <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:X:Z
"

est_cubic <- modsem_2smm(m_cubic, data = oneInt, verbose = FALSE)
testthat::expect_s3_class(est_cubic, "modsem_2smm")

co_cubic <- coef(est_cubic)
testthat::expect_true("Y~X:X:Z" %in% names(co_cubic))
testthat::expect_true("Y~X" %in% names(co_cubic))
testthat::expect_true("Y~Z" %in% names(co_cubic))
testthat::expect_equal(length(co_cubic), 3L)

# X and Z main effects should still be recovered reasonably
testthat::expect_true(co_cubic[["Y~X"]] > 0.3)
testthat::expect_true(co_cubic[["Y~Z"]] > 0.3)

# ===========================================================================
# Test 8: Mixed degree model (X:Z + X:X + X:X:Z)
# ===========================================================================
m_mixed <- "
  X =~ x1 + x2 + x3
  Z =~ z1 + z2 + z3
  Y =~ y1 + y2 + y3
  Y ~ X + Z + X:Z + X:X + X:X:Z
"

est_mixed <- modsem_2smm(m_mixed, data = oneInt, verbose = FALSE)
co_mixed <- coef(est_mixed)

testthat::expect_true("Y~X:Z" %in% names(co_mixed))
testthat::expect_true("Y~X:X" %in% names(co_mixed))
testthat::expect_true("Y~X:X:Z" %in% names(co_mixed))
testthat::expect_equal(length(co_mixed), 5L)

# X:Z should still be the dominant interaction effect
testthat::expect_true(abs(co_mixed[["Y~X:Z"]]) > abs(co_mixed[["Y~X:X"]]))
testthat::expect_true(abs(co_mixed[["Y~X:Z"]]) > abs(co_mixed[["Y~X:X:Z"]]))

# ===========================================================================
# Test 9: Cubic model with error.dist = "normal" works
# ===========================================================================
est_cubic_norm <- modsem_2smm(m_cubic, data = oneInt,
                               error.dist = "normal", verbose = FALSE)
testthat::expect_s3_class(est_cubic_norm, "modsem_2smm")
testthat::expect_true("Y~X:X:Z" %in% names(coef(est_cubic_norm)))

# ===========================================================================
# Test 10: Helper function - term_to_multiindex
# ===========================================================================
fnames <- c("X", "Z", "Y")

r1 <- term_to_multiindex("(Intercept)", fnames)
testthat::expect_equal(unname(r1), c(0L, 0L, 0L))

r2 <- term_to_multiindex("X", fnames)
testthat::expect_equal(r2[["X"]], 1L)
testthat::expect_equal(r2[["Z"]], 0L)

r3 <- term_to_multiindex("X:Z", fnames)
testthat::expect_equal(r3[["X"]], 1L)
testthat::expect_equal(r3[["Z"]], 1L)
testthat::expect_equal(r3[["Y"]], 0L)

r4 <- term_to_multiindex("X:X:Z", fnames)
testthat::expect_equal(r4[["X"]], 2L)
testthat::expect_equal(r4[["Z"]], 1L)

r5 <- term_to_multiindex("X:X:X", fnames)
testthat::expect_equal(r5[["X"]], 3L)
testthat::expect_equal(r5[["Z"]], 0L)

# ===========================================================================
# Test 11: Helper function - moment_from_cumulants (normal distribution)
# ===========================================================================
# For N(0, sigma^2): kappa_1=0, kappa_2=sigma^2, kappa_k=0 for k>=3
sigma2 <- 2
kappa_normal <- c(0, sigma2, rep(0, 8))

# mu_2 = sigma^2
testthat::expect_equal(moment_from_cumulants(2, kappa_normal), sigma2)

# mu_4 = 3*sigma^4
testthat::expect_equal(moment_from_cumulants(4, kappa_normal), 3 * sigma2^2)

# mu_6 = 15*sigma^6
testthat::expect_equal(moment_from_cumulants(6, kappa_normal), 15 * sigma2^3)

# Odd moments are zero
testthat::expect_equal(moment_from_cumulants(1, kappa_normal), 0)
testthat::expect_equal(moment_from_cumulants(3, kappa_normal), 0)
testthat::expect_equal(moment_from_cumulants(5, kappa_normal), 0)

# ===========================================================================
# Test 12: Helper function - set_partitions_no_singletons counts
# ===========================================================================
testthat::expect_equal(length(set_partitions_no_singletons(2)), 1L)
testthat::expect_equal(length(set_partitions_no_singletons(3)), 1L)
testthat::expect_equal(length(set_partitions_no_singletons(4)), 4L)
testthat::expect_equal(length(set_partitions_no_singletons(5)), 11L)
testthat::expect_equal(length(set_partitions_no_singletons(6)), 41L)

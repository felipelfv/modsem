#' Two-Stage Method of Moments (2SMM) for Latent Interactions
#'
#' @description
#' `modsem_2smm()` estimates latent variable interaction models using the
#' Two-Stage Method of Moments approach (Wall & Amemiya 2000, 2003).
#'
#' Stage 1 fits a CFA via lavaan and computes Bartlett factor scores.
#' Stage 2 performs bias-corrected method-of-moments regression on the
#' factor scores to estimate structural parameters including interaction effects.
#'
#' @param model.syntax lavaan-style model syntax
#' @param data A data frame with observed variables
#' @param center.data Center observed variables before estimation?
#' @param standardize.data Standardize observed variables before estimation?
#' @param se.type Type of standard errors: `"sandwich"` (default, HC0)
#' @param error.dist Distribution assumption for measurement errors:
#'   `"general"` (default, distribution-free) or `"normal"` (assumes normal errors)
#' @param verbose Print progress?
#' @param ... Additional arguments passed to `lavaan::cfa()` in Stage 1
#'
#' @return An object of class `c("modsem_2smm", "modsem")`
#'
#' @references
#' Wall, M. M. & Amemiya, Y. (2000). Estimation for polynomial structural
#' equation models. *Journal of the American Statistical Association*, 95(451),
#' 929-940.
#'
#' Wall, M. M. & Amemiya, Y. (2003). A method of moments technique for fitting
#' interaction effects in structural equation models. *British Journal of
#' Mathematical and Statistical Psychology*, 56(1), 47-63.
#'
#' @examples
#' \dontrun{
#' library(modsem)
#'
#' m1 <- "
#'   X =~ x1 + x2 + x3
#'   Z =~ z1 + z2 + z3
#'   Y =~ y1 + y2 + y3
#'   Y ~ X + Z + X:Z
#' "
#'
#' est <- modsem_2smm(m1, data = oneInt)
#' summary(est)
#' }
#'
#' @export
modsem_2smm <- function(model.syntax,
                         data,
                         center.data = FALSE,
                         standardize.data = FALSE,
                         se.type = "sandwich",
                         error.dist = "general",
                         verbose = TRUE,
                         ...) {

  if (is.null(model.syntax)) stop2("No model.syntax provided")
  if (!is.character(model.syntax) || length(model.syntax) != 1)
    stop2("model.syntax must be a single character string")
  if (is.null(data)) stop2("No data provided")

  data <- as.data.frame(data)
  error.dist <- match.arg(error.dist, c("general", "normal"))

  # Preprocess data
  if (center.data || standardize.data) {
    numeric_cols <- vapply(data, is.numeric, logical(1))
    if (center.data)
      data[numeric_cols] <- lapply(data[numeric_cols],
                                   function(x) x - mean(x, na.rm = TRUE))
    if (standardize.data)
      data[numeric_cols] <- lapply(data[numeric_cols], scale)
  }

  # Parse model syntax
  parTable <- modsemify(model.syntax)

  # Extract model components
  etas <- getSortedEtas(parTable, isLV = TRUE)
  xis  <- getXis(parTable, etas = etas, isLV = TRUE)
  lvs  <- getLVs(parTable)
  intTerms <- getIntTerms(parTable)

  if (length(intTerms) == 0)
    warning2("No interaction terms found in model syntax. ",
             "2SMM is designed for models with latent interactions.")

  if (length(etas) == 0) stop2("No endogenous latent variables found.")

  # Build CFA syntax (measurement model only)
  cfaParTable <- parTable[parTable$op == "=~", ]
  syntaxCFA   <- parTableToSyntax(cfaParTable)

  # Compute max polynomial degree from interaction terms
  if (length(intTerms) > 0) {
    max_degree <- max(vapply(intTerms, function(t) {
      length(strsplit(t, ":")[[1]])
    }, integer(1)))
  } else {
    max_degree <- 1L
  }

  if (verbose) {
    cat("modsem_2smm: Two-Stage Method of Moments\n")
    cat("==========================================\n")
    cat("  Exogenous LVs: ", paste(xis, collapse = ", "), "\n")
    cat("  Endogenous LVs:", paste(etas, collapse = ", "), "\n")
    cat("  Interactions:  ", paste(intTerms, collapse = ", "), "\n")
    cat("  Max poly degree:", max_degree, "\n")
    cat("  Error dist:    ", error.dist, "\n")
    cat("  N:             ", nrow(data), "\n\n")
  }

  # ---- Stage 1: CFA + Bartlett factor scores ----
  if (verbose) cat("Stage 1: Fitting CFA and computing Bartlett factor scores...\n")

  stage1 <- estimate2SMM_stage1(syntaxCFA, data = data, lvs = lvs, ...)

  # ---- Higher-order error moments ----
  if (verbose) cat("Estimating measurement error moments (", error.dist, ")...\n")

  max_order <- 2L * max_degree
  errMom <- estimateErrorMoments(
    W = stage1$W, Lambda = stage1$Lambda,
    Theta = stage1$Theta, v = stage1$v,
    factor_names = stage1$factor_names,
    max_order = max_order,
    error.dist = error.dist
  )

  # Create memoized error_moment_fn closure: given a multi-index s (named
  # integer vector), return E[prod e_l^{s_l}] using the estimated
  # eps_cumulants and W. Results are cached to avoid recomputing expensive
  # set-partition-based joint moments.
  W_mat <- stage1$W
  eps_cumulants <- errMom$eps_cumulants
  factor_names_ordered <- stage1$factor_names
  .emf_cache <- new.env(hash = TRUE, parent = emptyenv())

  error_moment_fn <- function(s) {
    key <- paste(s, collapse = ",")
    if (exists(key, envir = .emf_cache)) {
      return(get(key, envir = .emf_cache))
    }
    # Expand multi-index s into a vector of factor indices
    idx_vec <- integer(0)
    for (l in seq_along(s)) {
      if (s[l] > 0) {
        idx_vec <- c(idx_vec, rep(l, s[l]))
      }
    }
    val <- compute_error_joint_moment(idx_vec, W_mat, eps_cumulants)
    assign(key, val, envir = .emf_cache)
    val
  }

  # ---- Stage 2: Bias-corrected regression for each eta ----
  if (verbose) cat("Stage 2: Bias-corrected method-of-moments regression...\n")

  structExprs <- parTable[parTable$op == "~", ]
  all_alpha    <- list()
  all_se_info  <- list()

  for (eta in etas) {
    # Predictors for this eta
    eta_preds <- structExprs[structExprs$lhs == eta, "rhs"]

    linear_preds <- eta_preds[!grepl(":", eta_preds)]
    eta_ints     <- eta_preds[grepl(":", eta_preds)]

    if (length(linear_preds) == 0 && length(eta_ints) == 0) next

    stage2 <- estimate2SMM_stage2(
      f_hat = stage1$f_hat,
      Sigma_ee = stage1$Sigma_ee,
      eta_name = eta,
      linear_preds = linear_preds,
      int_terms = eta_ints,
      error_moment_fn = error_moment_fn,
      factor_names = factor_names_ordered
    )

    se_info <- compute2SMMSE(
      alpha_hat = stage2$alpha_hat,
      M_hat_inv = stage2$M_hat_inv,
      R = stage2$R,
      f_eta = stage2$f_eta,
      n = stage1$n
    )

    all_alpha[[eta]] <- list(
      alpha_hat = stage2$alpha_hat,
      linear_preds = linear_preds,
      int_terms = eta_ints,
      M_hat = stage2$M_hat,
      m_hat = stage2$m_hat,
      M_hat_inv = stage2$M_hat_inv,
      residuals = se_info$residuals
    )
    all_se_info[[eta]] <- se_info
  }

  # ---- Build parameter table ----
  if (verbose) cat("Building parameter table...\n")

  # Start with CFA parameter table
  cfa_pt <- lavaan::parameterEstimates(stage1$cfa_fit)
  cfa_pt <- cfa_pt[, c("lhs", "op", "rhs", "est", "se", "z", "pvalue",
                        "ci.lower", "ci.upper")]
  names(cfa_pt) <- c("lhs", "op", "rhs", "est", "std.error", "z.value",
                      "p.value", "ci.lower", "ci.upper")

  struct_rows <- NULL
  intercept_rows <- NULL
  resid_var_rows <- NULL

  for (eta in names(all_alpha)) {
    info <- all_alpha[[eta]]
    se_i <- all_se_info[[eta]]

    alpha <- info$alpha_hat
    se_alpha <- se_i$se

    # Structural regressions (exclude intercept)
    rhs_names <- c(info$linear_preds, info$int_terms)
    if (length(rhs_names) > 0) {
      sr <- data.frame(
        lhs = eta, op = "~", rhs = rhs_names,
        est = unname(alpha[-1]), std.error = unname(se_alpha[-1]),
        stringsAsFactors = FALSE
      )
      sr$z.value  <- sr$est / sr$std.error
      sr$p.value  <- 2 * stats::pnorm(-abs(sr$z.value))
      sr$ci.lower <- sr$est - CI_WIDTH * sr$std.error
      sr$ci.upper <- sr$est + CI_WIDTH * sr$std.error
      struct_rows <- rbind(struct_rows, sr)
    }

    # Intercept
    ir <- data.frame(
      lhs = eta, op = "~1", rhs = "",
      est = unname(alpha["(Intercept)"]),
      std.error = unname(se_alpha["(Intercept)"]),
      stringsAsFactors = FALSE
    )
    ir$z.value  <- ir$est / ir$std.error
    ir$p.value  <- 2 * stats::pnorm(-abs(ir$z.value))
    ir$ci.lower <- ir$est - CI_WIDTH * ir$std.error
    ir$ci.upper <- ir$est + CI_WIDTH * ir$std.error
    intercept_rows <- rbind(intercept_rows, ir)

    # Residual variance (bias-corrected)
    raw_resid_var <- mean(info$residuals^2)
    psi_hat <- max(raw_resid_var - stage1$Sigma_ee[eta, eta], 0)

    rv <- data.frame(
      lhs = eta, op = "~~", rhs = eta,
      est = psi_hat, std.error = NA_real_,
      z.value = NA_real_, p.value = NA_real_,
      ci.lower = NA_real_, ci.upper = NA_real_,
      stringsAsFactors = FALSE
    )
    resid_var_rows <- rbind(resid_var_rows, rv)
  }

  parTable_out <- rbind(cfa_pt, struct_rows, intercept_rows, resid_var_rows)
  rownames(parTable_out) <- NULL
  parTable_out <- fillColsParTable(parTable_out)
  parTable_out$label[is.na(parTable_out$label)] <- ""

  # Build coef vector and vcov for structural parameters
  coef_struct <- unlist(lapply(all_alpha, function(x) x$alpha_hat))
  vcov_struct <- as.matrix(Matrix::bdiag(lapply(all_se_info, function(x) x$vcov)))
  rownames(vcov_struct) <- colnames(vcov_struct) <- names(coef_struct)

  if (verbose) cat("\nEstimation complete. Use summary() for results.\n")

  # ---- Build return object ----
  out <- list(
    parTable         = parTable_out,
    coef.structural  = coef_struct,
    vcov.structural  = vcov_struct,
    cfa              = stage1$cfa_fit,
    factor.scores    = stage1$f_hat,
    Sigma_ee         = stage1$Sigma_ee,
    syntax           = model.syntax,
    data             = data,
    N                = stage1$n,
    method           = "2smm",
    se.type          = se.type,
    error.dist       = error.dist,
    etas             = etas,
    xis              = xis,
    intTerms         = intTerms,
    stage2           = all_alpha
  )

  class(out) <- c("modsem_2smm", "modsem")
  out
}


# =============================================================================
# S3 methods
# =============================================================================

#' @export
summary.modsem_2smm <- function(object, digits = 3, scientific = FALSE,
                                 ci = FALSE, ...) {
  out <- list(
    object = object,
    digits = digits,
    scientific = scientific,
    ci = ci
  )
  class(out) <- "summary_modsem_2smm"
  out
}


#' @export
print.summary_modsem_2smm <- function(x, ...) {
  object     <- x$object
  digits     <- x$digits
  scientific <- x$scientific
  ci         <- x$ci

  cat("\n")
  cat("modsem_2smm: Two-Stage Method of Moments\n")
  cat("==========================================\n")
  cat("  Method:     2SMM (Wall & Amemiya 2000, 2003)\n")
  cat("  Error dist:", object$error.dist, "\n")
  cat("  SE type:   ", object$se.type, "\n")
  cat("  N:         ", object$N, "\n")

  # CFA fit info
  cfa_fit <- object$cfa
  fit_measures <- lavaan::fitMeasures(cfa_fit, c("chisq", "df", "pvalue",
                                                  "cfi", "rmsea", "srmr"))
  cat("\nCFA Fit (Stage 1):\n")
  cat(sprintf("  Chi-sq(%d) = %.3f, p = %.3f\n",
              fit_measures["df"], fit_measures["chisq"], fit_measures["pvalue"]))
  cat(sprintf("  CFI = %.3f, RMSEA = %.3f, SRMR = %.3f\n",
              fit_measures["cfi"], fit_measures["rmsea"], fit_measures["srmr"]))

  cat("\nParameter Estimates (Stage 2):\n")

  printParTable(object$parTable, digits = digits, scientific = scientific,
                ci = ci, loadings = TRUE, regressions = TRUE,
                covariances = TRUE, intercepts = TRUE, variances = TRUE)

  invisible(x)
}


#' @export
print.modsem_2smm <- function(x, ...) {
  cat("modsem_2smm object\n")
  cat("Method: Two-Stage Method of Moments (Wall & Amemiya 2000, 2003)\n")
  cat("N:", x$N, "\n")
  cat("Error dist:", x$error.dist, "\n")
  cat("\nUse summary() for detailed results.\n")
  invisible(x)
}


#' @export
#' @importFrom stats coef
coef.modsem_2smm <- function(object, ...) {
  pt <- object$parTable
  struct <- pt[pt$op == "~", ]
  est <- struct$est
  names(est) <- paste0(struct$lhs, "~", struct$rhs)
  est
}


#' @export
#' @importFrom stats vcov
vcov.modsem_2smm <- function(object, ...) {
  object$vcov.structural
}


#' @export
#' @importFrom stats nobs
nobs.modsem_2smm <- function(object, ...) {
  object$N
}


#' @describeIn parameter_estimates Get parameter estimates of a
#'   \code{modsem_2smm} object
#' @export
parameter_estimates.modsem_2smm <- function(object, ...) {
  modsemParTable(object$parTable)
}


#' @describeIn standardized_estimates Method for \code{modsem_2smm} objects
#' @export
standardized_estimates.modsem_2smm <- function(object, ...) {
  pt <- parameter_estimates(object)
  standardized_estimates.data.frame(pt, ...)
}


#' @export
var_interactions.modsem_2smm <- function(object, ...) {
  pt <- parameter_estimates(object)
  var_interactions.data.frame(pt, ...)
}


#' @export
is_interaction_model.modsem_2smm <- function(object) {
  length(object$intTerms) > 0
}


#' @export
modsem_inspect.modsem_2smm <- function(object, what = NULL, ...) {
  if (is.null(what)) what <- "default"

  switch(what,
    "default"       = , "partable" = object$parTable,
    "cfa"           = object$cfa,
    "sigma_ee"      = object$Sigma_ee,
    "factor.scores" = object$factor.scores,
    "vcov"          = object$vcov.structural,
    "coef"          = coef(object),
    "N"             = , "nobs" = object$N,
    "method"        = object$method,
    stop2("Unknown 'what': ", what)
  )
}

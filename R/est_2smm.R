# ===========================================================================
# Core estimation functions for the 2SMM (Two-Stage Method of Moments)
# approach (Wall & Amemiya 2000, 2003), generalized to arbitrary polynomial
# degree (cubic, quartic, etc.)
# ===========================================================================


# Stage 1: CFA estimation + Bartlett factor scores
estimate2SMM_stage1 <- function(syntaxCFA, data, lvs, ...) {
  cfa_fit <- lavaan::cfa(syntaxCFA, data = data, meanstructure = TRUE, ...)

  if (!lavaan::lavInspect(cfa_fit, "converged")) {
    stop2("CFA did not converge. Cannot compute Bartlett factor scores.")
  }

  est <- lavaan::lavInspect(cfa_fit, "est")
  Lambda <- est$lambda     # p x q (indicators x factors)
  tau <- as.vector(est$nu) # p-vector (indicator intercepts)
  Theta <- est$theta       # p x p (indicator error covariance)

  factor_names <- colnames(Lambda)
  indicator_names <- rownames(Lambda)
  p <- nrow(Lambda)
  q <- ncol(Lambda)

  # Check for non-diagonal Theta
  off_diag <- Theta[lower.tri(Theta)]
  if (any(abs(off_diag) > 1e-10)) {
    warning2("Theta has off-diagonal elements. ",
             "Higher-order moment estimation assumes independent ",
             "measurement errors.")
  }

  # Bartlett weight matrix: W = (Lambda'Theta^{-1}Lambda)^{-1} Lambda'Theta^{-1}
  Theta_inv <- solve(Theta)
  LtTiL <- crossprod(Lambda, Theta_inv %*% Lambda)
  Sigma_ee <- solve(LtTiL)
  W <- Sigma_ee %*% crossprod(Lambda, Theta_inv)

  # Factor scores: f_hat_t = W(z_t - tau)
  Z <- as.matrix(data[, indicator_names, drop = FALSE])
  Z_centered <- sweep(Z, 2, tau)
  f_hat <- Z_centered %*% t(W)
  colnames(f_hat) <- factor_names

  # CFA residuals: v_t = z_t - tau - Lambda * f_hat_t
  v <- Z_centered - f_hat %*% t(Lambda)

  # Reorder to match lvs ordering
  f_hat <- f_hat[, lvs, drop = FALSE]
  Sigma_ee <- Sigma_ee[lvs, lvs, drop = FALSE]
  W <- W[lvs, , drop = FALSE]

  list(
    cfa_fit = cfa_fit,
    Lambda = Lambda,
    tau = tau,
    Theta = Theta,
    Theta_inv = Theta_inv,
    f_hat = f_hat,
    Sigma_ee = Sigma_ee,
    W = W,
    v = v,
    factor_names = lvs,
    indicator_names = indicator_names,
    n = nrow(Z)
  )
}


# ===========================================================================
# Mathematical helper functions for generalized moment computation
# ===========================================================================

# Convert a term string to a multi-index vector over factor_names
# e.g., "X:X:Z" with factor_names = c("X","Z","Y") -> c(X=2, Z=1, Y=0)
#        "(Intercept)" -> c(X=0, Z=0, Y=0)
#        "X" -> c(X=1, Z=0, Y=0)
term_to_multiindex <- function(term, factor_names) {
  r <- integer(length(factor_names))
  names(r) <- factor_names
  if (term == "(Intercept)") return(r)
  vars <- strsplit(term, ":")[[1]]
  for (v in vars) {
    r[v] <- r[v] + 1L
  }
  r
}


# Univariate moment-cumulant recursion:
# mu_n = sum_{j=1}^{n} C(n-1,j-1) * kappa_j * mu_{n-j}
# where mu_0 = 1, kappa_vec[k] = kappa_k
moment_from_cumulants <- function(k, kappa_vec) {
  if (k == 0) return(1)
  mu <- numeric(k + 1)
  mu[1] <- 1  # mu_0
  for (n in seq_len(k)) {
    s <- 0
    for (j in seq_len(n)) {
      if (j <= length(kappa_vec)) {
        s <- s + choose(n - 1, j - 1) * kappa_vec[j] * mu[n - j + 1]
      }
    }
    mu[n + 1] <- s
  }
  mu[k + 1]
}


# Module-level cache for set partitions
.partition_cache <- new.env(hash = TRUE, parent = emptyenv())

# Enumerate all set partitions of {1,...,n} where every block has size >= 2
# Returns a list of partitions, each partition is a list of integer vectors (blocks)
set_partitions_no_singletons <- function(n) {
  key <- as.character(n)
  if (exists(key, envir = .partition_cache)) {
    return(get(key, envir = .partition_cache))
  }

  result <- .enumerate_partitions(seq_len(n))
  assign(key, result, envir = .partition_cache)
  result
}


# Recursive helper: enumerate partitions of 'elements' with no singletons
.enumerate_partitions <- function(elements) {
  n <- length(elements)
  if (n == 0) return(list(list()))  # single empty partition
  if (n == 1) return(list())        # impossible: one element can't form block >=2

  # Take first element; it must join with at least one other
  first <- elements[1]
  rest <- elements[-1]

  result <- list()

  # Choose which other elements join first in a block: block size from 2..n
  for (block_size in seq(2, n, by = 1)) {
    remainder <- n - block_size
    if (remainder == 1) next  # remainder of 1 can't form block >=2

    # All combinations of (block_size-1) elements from rest
    combos <- utils::combn(rest, block_size - 1, simplify = FALSE)
    for (combo in combos) {
      block <- c(first, combo)
      remaining <- setdiff(rest, combo)
      sub_parts <- .enumerate_partitions(remaining)
      for (sp in sub_parts) {
        result <- c(result, list(c(list(block), sp)))
      }
    }
  }

  result
}


# Compute the joint moment E[prod_l e_{idx[l]}] for measurement error e = W * eps
# using the moment-cumulant formula (set partitions with no singletons since
# kappa_1(e) = 0 for Bartlett scores centered at 0)
#
# idx_vec: integer vector of factor indices (e.g., c(1,1,3) for e_1*e_1*e_3)
# W: q x p weight matrix
# eps_cumulants: list where eps_cumulants[[k]] is a p-vector of k-th cumulants (k>=2)
compute_error_joint_moment <- function(idx_vec, W, eps_cumulants) {
  n <- length(idx_vec)
  if (n == 0) return(1)
  if (n == 1) return(0)  # E[e_a] = 0

  partitions <- set_partitions_no_singletons(n)
  if (length(partitions) == 0) return(0)

  p <- ncol(W)
  total <- 0

  for (partition in partitions) {
    # Product over blocks
    block_prod <- 1
    for (block in partition) {
      k <- length(block)
      if (k > length(eps_cumulants) || is.null(eps_cumulants[[k]])) {
        # Cumulant not available; treat as 0
        block_prod <- 0
        break
      }
      kappa_k <- eps_cumulants[[k]]  # p-vector of k-th cumulants
      # Joint cumulant of (e_{idx[b1]}, ..., e_{idx[bk]}) =
      # sum_j prod_{l in block} W[idx[l], j] * kappa_k[j]
      factor_indices <- idx_vec[block]
      # Compute W-product for each indicator j
      w_prod <- rep(1, p)
      for (fi in factor_indices) {
        w_prod <- w_prod * W[fi, ]
      }
      block_prod <- block_prod * sum(w_prod * kappa_k)
    }
    total <- total + block_prod
  }

  total
}


# Compute raw sample moment: mean(prod_{l: r_l > 0} f_hat[,l]^r_l)
compute_raw_moment <- function(f_hat, r) {
  active <- which(r > 0)
  if (length(active) == 0) return(1)
  n <- nrow(f_hat)
  prod_vec <- rep(1, n)
  for (l in active) {
    prod_vec <- prod_vec * f_hat[, l]^r[l]
  }
  mean(prod_vec)
}


# Compute unbiased moment U_r(f_hat) using recursive bias correction:
# U_r = raw_moment(f_hat, r) - SUM_{0 < s <= r} C(r,s) * mu_e(s) * U_{r-s}
# where mu_e(s) is the joint moment of measurement errors e for multi-index s
#
# f_hat: n x q matrix of factor scores
# r: named integer multi-index
# error_moment_fn: function(s) returning E[prod e_l^{s_l}]
# memo: environment for memoization
unbiased_moment <- function(r, f_hat, error_moment_fn, memo) {
  key <- paste(r, collapse = ",")
  if (exists(key, envir = memo)) return(get(key, envir = memo))

  total_order <- sum(r)
  if (total_order == 0) {
    assign(key, 1, envir = memo)
    return(1)
  }

  raw <- compute_raw_moment(f_hat, r)

  # Enumerate all sub-indices s with 0 <= s_l <= r_l, excluding s = 0
  q <- length(r)
  grid_list <- lapply(r, function(ri) 0:ri)
  all_s <- as.matrix(expand.grid(grid_list))
  colnames(all_s) <- names(r)

  correction <- 0
  for (row_idx in seq_len(nrow(all_s))) {
    s <- all_s[row_idx, ]
    s_sum <- sum(s)
    if (s_sum == 0) next          # skip s = 0
    if (s_sum == 1) next          # E[e_a] = 0 always
    if (all(s == r)) {
      # When s == r, U_{r-s} = U_0 = 1
      mu_e_s <- error_moment_fn(s)
      if (abs(mu_e_s) < 1e-30) next
      # C(r, s) = prod C(r_l, s_l) = 1 when s == r
      correction <- correction + mu_e_s
      next
    }

    mu_e_s <- error_moment_fn(s)
    if (abs(mu_e_s) < 1e-30) next

    # Multinomial coefficient C(r, s) = prod_l C(r_l, s_l)
    C_rs <- prod(vapply(seq_len(q), function(l) choose(r[l], s[l]), numeric(1)))

    r_minus_s <- r - s
    U_rms <- unbiased_moment(r_minus_s, f_hat, error_moment_fn, memo)

    correction <- correction + C_rs * mu_e_s * U_rms
  }

  result <- raw - correction
  assign(key, result, envir = memo)
  result
}


# ===========================================================================
# Generalized error moment estimation
# ===========================================================================

# Estimate cumulants of epsilon (indicator errors) up to max_order
# Returns eps_cumulants: list where eps_cumulants[[k]] is a p-vector (k = 2, 3, ...)
estimateErrorMoments <- function(W, Lambda, Theta, v, factor_names,
                                 max_order = 4, error.dist = "general") {
  q <- nrow(W)
  p <- ncol(W)
  n <- nrow(v)
  theta_diag <- diag(Theta)

  # kappa_2 is always known from CFA: Var(eps_j) = Theta_jj
  eps_cumulants <- list()
  eps_cumulants[[2]] <- theta_diag

  if (error.dist == "normal" || max_order < 3) {
    # All cumulants >= 3 are zero for normal distribution
    if (max_order >= 3) {
      for (k in 3:max_order) {
        eps_cumulants[[k]] <- rep(0, p)
      }
    }
    return(list(eps_cumulants = eps_cumulants))
  }

  # General case: estimate from CFA residuals
  # v = C * epsilon where C = I_p - Lambda %*% W
  C_mat <- diag(p) - Lambda %*% W

  # Iteratively estimate kappa_k for k = 3, 4, ..., max_order
  # At each step, compute E[v_i^k] and subtract known lower-order contributions
  for (k in 3:max_order) {
    s_k <- colMeans(v^k)

    # Compute known contribution from lower-order cumulants
    # E[v_i^k] = E[(sum_j C[i,j] * eps_j)^k]
    # This is a sum over set partitions of {1,...,k}:
    # For each partition pi, product over blocks B:
    #   sum_j C[i,j]^|B| * kappa_{|B|}(eps_j)
    # We need to separate the contribution of kappa_k (unknown) from known lower terms
    known_lower <- numeric(p)
    for (i in seq_len(p)) {
      known_lower[i] <- .compute_univariate_mixed_moment(
        k, C_mat[i, ], eps_cumulants, exclude_order = k
      )
    }

    # Design: contribution of kappa_k to E[v_i^k]
    # Only partition with single block {1,...,k}: sum_j C[i,j]^k * kappa_k[j]
    A_k <- C_mat^k

    s_k_adj <- s_k - known_lower

    kappa_k <- tryCatch(
      as.vector(solve(crossprod(A_k), crossprod(A_k, s_k_adj))),
      error = function(e) {
        warning2("Could not estimate order-", k, " cumulants of measurement errors. ",
                 "Falling back to 0.")
        rep(0, p)
      }
    )

    eps_cumulants[[k]] <- kappa_k
  }

  list(eps_cumulants = eps_cumulants)
}


# Compute E[(sum_j c_j * eps_j)^k] using moment-cumulant formula,
# where eps_j are independent with known cumulants, optionally excluding
# one cumulant order (set to 0) for the purpose of solving for it
.compute_univariate_mixed_moment <- function(k, c_vec, eps_cumulants,
                                              exclude_order = NULL) {
  # For a linear combination Y = sum c_j eps_j of independent variables,
  # kappa_m(Y) = sum_j c_j^m * kappa_m(eps_j)
  # Then moment of Y from its cumulants via moment-cumulant recursion

  p <- length(c_vec)

  # Compute cumulants of Y up to order k
  y_cumulants <- numeric(k)
  for (m in seq_len(k)) {
    if (m <= length(eps_cumulants) && !is.null(eps_cumulants[[m]])) {
      if (!is.null(exclude_order) && m == exclude_order) {
        y_cumulants[m] <- 0
      } else {
        y_cumulants[m] <- sum(c_vec^m * eps_cumulants[[m]])
      }
    }
  }

  # Moment from cumulants
  moment_from_cumulants(k, y_cumulants)
}


# ===========================================================================
# Stage 2: Generalized bias-corrected method-of-moments regression
# ===========================================================================

estimate2SMM_stage2 <- function(f_hat, Sigma_ee, eta_name,
                                linear_preds, int_terms,
                                error_moment_fn, factor_names) {
  n <- nrow(f_hat)
  q <- length(factor_names)

  f_eta <- f_hat[, eta_name]

  # Build regressor names and multi-indices
  reg_names <- c("(Intercept)", linear_preds, int_terms)
  n_reg <- length(reg_names)
  n_lin <- length(linear_preds)
  n_int <- length(int_terms)

  # Multi-index for each regressor
  reg_multiidx <- lapply(reg_names, term_to_multiindex, factor_names = factor_names)
  names(reg_multiidx) <- reg_names

  # Multi-index for eta
  eta_midx <- integer(q)
  names(eta_midx) <- factor_names
  eta_midx[eta_name] <- 1L

  # Build regressor matrix R = [1, f_xi1, ..., f_xia*f_xib*..., ...]
  R <- matrix(0, nrow = n, ncol = n_reg)
  R[, 1] <- 1
  for (j in seq_len(n_reg)[-1]) {
    r_j <- reg_multiidx[[j]]
    active <- which(r_j > 0)
    if (length(active) == 0) {
      R[, j] <- 1
    } else {
      col <- rep(1, n)
      for (l in active) {
        col <- col * f_hat[, l]^r_j[l]
      }
      R[, j] <- col
    }
  }

  # Memoization environment for unbiased moments
  memo <- new.env(hash = TRUE, parent = emptyenv())

  # Compute M_hat[i,j] = U_{r_i + r_j}
  M_hat <- matrix(0, nrow = n_reg, ncol = n_reg)
  for (i in seq_len(n_reg)) {
    for (j in i:n_reg) {
      r_ij <- reg_multiidx[[i]] + reg_multiidx[[j]]
      val <- unbiased_moment(r_ij, f_hat, error_moment_fn, memo)
      M_hat[i, j] <- val
      M_hat[j, i] <- val
    }
  }

  # Compute m_hat[i] = U_{r_i + e_eta}
  m_hat <- numeric(n_reg)
  for (i in seq_len(n_reg)) {
    r_i_eta <- reg_multiidx[[i]] + eta_midx
    m_hat[i] <- unbiased_moment(r_i_eta, f_hat, error_moment_fn, memo)
  }

  # Solve: alpha_hat = M_hat^{-1} m_hat
  M_hat_inv <- tryCatch(
    solve(M_hat),
    error = function(e) {
      warning2("M_hat is singular or near-singular. Using generalized inverse.")
      MASS::ginv(M_hat)
    }
  )

  alpha_hat <- as.vector(M_hat_inv %*% m_hat)
  names(alpha_hat) <- reg_names

  list(
    alpha_hat = alpha_hat,
    reg_names = reg_names,
    M_hat = M_hat,
    m_hat = m_hat,
    M_hat_inv = M_hat_inv,
    R = R,
    f_eta = f_eta,
    n_lin = n_lin,
    n_int = n_int
  )
}


# ===========================================================================
# Standard errors & parameter table (unchanged)
# ===========================================================================

# Compute sandwich standard errors for 2SMM
compute2SMMSE <- function(alpha_hat, M_hat_inv, R, f_eta, n) {
  # Residuals
  e_hat <- f_eta - R %*% alpha_hat

  # Observation-level estimating equation: ell_t = e_t * r_t
  # Omega* = (1/n) sum(e_t^2 * r_t r_t')
  ell <- as.vector(e_hat) * R
  Omega_star <- crossprod(ell) / n

  # V(alpha) = (1/n) M_hat^{-1} Omega* M_hat^{-1}
  vcov_alpha <- (M_hat_inv %*% Omega_star %*% M_hat_inv) / n
  se <- sqrt(pmax(diag(vcov_alpha), 0))
  names(se) <- names(alpha_hat)

  list(vcov = vcov_alpha, se = se, Omega_star = Omega_star,
       residuals = as.vector(e_hat))
}


# Build parameter table combining CFA and structural estimates
build2SMMParTable <- function(cfa_fit, alpha_hat, se_alpha, eta_name,
                              linear_preds, int_terms, Sigma_ee,
                              residuals_structural) {
  # CFA parameter estimates from lavaan
  cfa_pt <- lavaan::parameterEstimates(cfa_fit)

  # Keep measurement model + indicator variances/covariances + factor covariances
  cfa_pt <- cfa_pt[, c("lhs", "op", "rhs", "est", "se", "z", "pvalue",
                        "ci.lower", "ci.upper")]
  names(cfa_pt) <- c("lhs", "op", "rhs", "est", "std.error", "z.value",
                      "p.value", "ci.lower", "ci.upper")

  # Structural parameters
  n <- length(residuals_structural)

  struct_rows <- data.frame(
    lhs = eta_name,
    op = "~",
    rhs = c(linear_preds, int_terms),
    est = alpha_hat[-1],  # exclude intercept
    std.error = se_alpha[-1],
    stringsAsFactors = FALSE
  )
  struct_rows$z.value <- struct_rows$est / struct_rows$std.error
  struct_rows$p.value <- 2 * stats::pnorm(-abs(struct_rows$z.value))
  struct_rows$ci.lower <- struct_rows$est - CI_WIDTH * struct_rows$std.error
  struct_rows$ci.upper <- struct_rows$est + CI_WIDTH * struct_rows$std.error

  # Intercept row
  intercept_row <- data.frame(
    lhs = eta_name,
    op = "~1",
    rhs = "",
    est = alpha_hat["(Intercept)"],
    std.error = se_alpha["(Intercept)"],
    stringsAsFactors = FALSE
  )
  intercept_row$z.value <- intercept_row$est / intercept_row$std.error
  intercept_row$p.value <- 2 * stats::pnorm(-abs(intercept_row$z.value))
  intercept_row$ci.lower <- intercept_row$est - CI_WIDTH * intercept_row$std.error
  intercept_row$ci.upper <- intercept_row$est + CI_WIDTH * intercept_row$std.error

  # Residual variance: psi = (1/n) sum(resid^2) - Sigma_ee[eta, eta]
  raw_resid_var <- mean(residuals_structural^2)
  psi_hat <- raw_resid_var - Sigma_ee[eta_name, eta_name]
  psi_hat <- max(psi_hat, 0)  # ensure non-negative

  resid_var_row <- data.frame(
    lhs = eta_name,
    op = "~~",
    rhs = eta_name,
    est = psi_hat,
    std.error = NA_real_,
    z.value = NA_real_,
    p.value = NA_real_,
    ci.lower = NA_real_,
    ci.upper = NA_real_,
    stringsAsFactors = FALSE
  )

  # Combine
  parTable <- rbind(cfa_pt, struct_rows, intercept_row, resid_var_row)
  rownames(parTable) <- NULL

  parTable
}

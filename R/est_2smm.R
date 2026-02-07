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
    Z = Z,
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
# using the recursive multivariate moment-cumulant formula.
#
# Since kappa_1(e) = 0 for Bartlett scores (centered), we use:
#   M(S) = SUM over subsets B of S\{first} with |B| >= 1:
#            C(|S|-2, |B|-1) * cum({first} ∪ B) * M(S \ B \ {first})
#            [then add the all-in-one-block term: cum(S)]
#
# More precisely, this is the standard recursion:
#   mu(S) = SUM_{B: first in B, |B| >= 2} C(|S|-2, |B|-2) * cum(B) * mu(S \ B)
#
# This avoids explicit set-partition enumeration (which is super-exponential).
# With memoization on subsets this is O(3^n) worst case, much better than
# enumerating all partitions for n >= 10.
#
# idx_vec: integer vector of factor indices (e.g., c(1,1,3) for e_1*e_1*e_3)
# W: q x p weight matrix
# eps_cumulants: list where eps_cumulants[[k]] is a p-vector of k-th cumulants (k>=2)
compute_error_joint_moment <- function(idx_vec, W, eps_cumulants) {
  n <- length(idx_vec)
  if (n == 0) return(1)
  if (n == 1) return(0)  # E[e_a] = 0

  p <- ncol(W)

  # Precompute W rows for each factor index used
  # W_rows[[i]] is the i-th row of W (a p-vector)
  # We use integer positions in idx_vec as keys
  unique_fi <- unique(idx_vec)
  W_rows <- list()
  for (fi in unique_fi) {
    W_rows[[as.character(fi)]] <- W[fi, ]
  }

  # Helper: compute joint cumulant cum(e_{idx[positions]})
  # = sum_j prod_{l in positions} W[idx[l], j] * kappa_{|positions|}[j]
  joint_cum <- function(positions) {
    k <- length(positions)
    if (k < 2) return(0)  # kappa_1 = 0
    if (k > length(eps_cumulants) || is.null(eps_cumulants[[k]])) return(0)
    kappa_k <- eps_cumulants[[k]]
    w_prod <- rep(1, p)
    for (pos in positions) {
      w_prod <- w_prod * W_rows[[as.character(idx_vec[pos])]]
    }
    sum(w_prod * kappa_k)
  }

  # Memoized recursive moment computation
  # Elements is an integer vector of positions into idx_vec
  memo_env <- new.env(hash = TRUE, parent = emptyenv())

  rec_moment <- function(elems) {
    ne <- length(elems)
    if (ne == 0) return(1)
    if (ne == 1) return(0)  # kappa_1 = 0

    key <- paste(elems, collapse = ",")
    if (exists(key, envir = memo_env)) return(get(key, envir = memo_env))

    first <- elems[1]
    rest <- elems[-1]
    nr <- length(rest)

    total <- 0

    # Enumerate subsets B_rest of rest with |B_rest| >= 1
    # Block = {first} union B_rest, so block size >= 2
    for (bsize_rest in seq_len(nr)) {
      # Remainder after removing B_rest from rest must be 0 or >= 2
      remainder_size <- nr - bsize_rest
      if (remainder_size == 1) next  # can't form partition with singleton remainder

      combos <- utils::combn(nr, bsize_rest, simplify = FALSE)
      for (combo_idx in combos) {
        block_positions <- c(first, rest[combo_idx])
        remaining <- rest[-combo_idx]

        # Joint cumulant of the block
        jc <- joint_cum(block_positions)
        if (abs(jc) < 1e-30) next

        # Combinatorial coefficient: C(ne-2, bsize_rest-1)
        # This counts the number of ways to reach this partition ordering
        coeff <- choose(ne - 2L, bsize_rest - 1L)

        # Recursive moment of the remaining elements
        m_rem <- rec_moment(remaining)

        total <- total + coeff * jc * m_rem
      }
    }

    assign(key, total, envir = memo_env)
    total
  }

  rec_moment(seq_len(n))
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


# Individual-level bias-corrected contributions (n-vector per multi-index r).
# Same recursion as unbiased_moment() but returns one value per observation
# instead of a scalar mean:
#   bc_t(r) = raw_product_t(r) - SUM_{0<s<=r, |s|>=2} C(r,s) * mu_e(s) * bc_t(r-s)
# For sum(r)=0: rep(1,n). For sum(r)=1: raw product (no correction).
compute_bc_individual <- function(r, f_hat, error_moment_fn, memo) {
  key <- paste(r, collapse = ",")
  if (exists(key, envir = memo)) return(get(key, envir = memo))

  n <- nrow(f_hat)
  total_order <- sum(r)

  if (total_order == 0) {
    result <- rep(1, n)
    assign(key, result, envir = memo)
    return(result)
  }

  # Raw product vector: prod_{l: r_l > 0} f_hat[,l]^r_l
  active <- which(r > 0)
  raw_vec <- rep(1, n)
  for (l in active) {
    raw_vec <- raw_vec * f_hat[, l]^r[l]
  }

  if (total_order == 1) {
    # No correction: E[e] = 0, so all corrections require |s| >= 2 > total_order
    assign(key, raw_vec, envir = memo)
    return(raw_vec)
  }

  # Enumerate all sub-indices s with 0 <= s_l <= r_l, excluding s = 0
  q <- length(r)
  grid_list <- lapply(r, function(ri) 0:ri)
  all_s <- as.matrix(expand.grid(grid_list))
  colnames(all_s) <- names(r)

  correction <- rep(0, n)
  for (row_idx in seq_len(nrow(all_s))) {
    s <- all_s[row_idx, ]
    s_sum <- sum(s)
    if (s_sum == 0) next
    if (s_sum == 1) next  # E[e_a] = 0

    mu_e_s <- error_moment_fn(s)
    if (abs(mu_e_s) < 1e-30) next

    if (all(s == r)) {
      # When s == r, bc_{r-s} = bc_0 = rep(1,n), C(r,s) = 1
      correction <- correction + mu_e_s
      next
    }

    C_rs <- prod(vapply(seq_len(q), function(l) choose(r[l], s[l]), numeric(1)))
    r_minus_s <- r - s
    bc_rms <- compute_bc_individual(r_minus_s, f_hat, error_moment_fn, memo)
    correction <- correction + C_rs * mu_e_s * bc_rms
  }

  result <- raw_vec - correction
  assign(key, result, envir = memo)
  result
}


# Compute bias-corrected ell matrix (Paper Eq 13-14):
#   ell_t[j] = m_hat_t[j] - sum_k M_hat_t[j,k] * alpha_hat[k]
# where m_hat_t[j] = bc_t(r_j + r_eta) and M_hat_t[j,k] = bc_t(r_j + r_k)
# are individual-level bias-corrected contributions.
# Returns n x n_reg matrix.
compute_bc_ell <- function(f_hat, alpha_hat, reg_multiidx,
                            eta_name, factor_names, error_moment_fn) {
  n     <- nrow(f_hat)
  q     <- length(factor_names)
  n_reg <- length(reg_multiidx)

  eta_midx <- integer(q)
  names(eta_midx) <- factor_names
  eta_midx[eta_name] <- 1L

  memo <- new.env(hash = TRUE, parent = emptyenv())

  ell <- matrix(0, nrow = n, ncol = n_reg)

  for (j in seq_len(n_reg)) {
    r_j <- reg_multiidx[[j]]

    # m_hat_t[j] = bc_t(r_j + r_eta)
    m_hat_j <- compute_bc_individual(r_j + eta_midx, f_hat, error_moment_fn, memo)

    # sum_k alpha_hat[k] * bc_t(r_j + r_k)
    M_alpha_j <- rep(0, n)
    for (k in seq_len(n_reg)) {
      bc_jk <- compute_bc_individual(r_j + reg_multiidx[[k]], f_hat,
                                      error_moment_fn, memo)
      M_alpha_j <- M_alpha_j + alpha_hat[k] * bc_jk
    }

    ell[, j] <- m_hat_j - M_alpha_j
  }

  ell
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

  # Proactive check: C has rank p - q, so C^k may be rank-deficient.
  # OLS for kappa_k solves a p x p system with design A_k = C^k.
  # When p - q < q (i.e., fewer than ~3 indicators per factor on average),
  # the system is likely underdetermined and kappa_k estimation will fail.
  if (p <= 2 * q) {
    warning2("Too few indicators for non-normal cumulant estimation ",
             "(", p, " indicators, ", q, " factors; need p >= 2*q = ", 2 * q, "). ",
             "Higher-order cumulants (kappa_3, kappa_4, ...) will fall back to 0 ",
             "(assumes normality). Consider using more indicators per factor ",
             "or error.dist = \"normal\".")
  }

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
        warning2("Could not estimate order-", k,
                 " cumulants of measurement errors ",
                 "(rank(C^k) too low; need more indicators per factor). ",
                 "Falling back to 0 (assumes normality for this order).")
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


# Factory function to create a memoized error moment function.
# Given W (weight matrix) and eps_cumulants (list of cumulant vectors),
# returns a function(s) that computes E[prod e_l^{s_l}] with caching.
make_error_moment_fn <- function(W, eps_cumulants) {
  cache <- new.env(hash = TRUE, parent = emptyenv())
  function(s) {
    key <- paste(s, collapse = ",")
    if (exists(key, envir = cache)) {
      return(get(key, envir = cache))
    }
    idx_vec <- integer(0)
    for (l in seq_along(s)) {
      if (s[l] > 0) {
        idx_vec <- c(idx_vec, rep(l, s[l]))
      }
    }
    val <- compute_error_joint_moment(idx_vec, W, eps_cumulants)
    assign(key, val, envir = cache)
    val
  }
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
# If C_hat and T_hat are provided, use the full sandwich:
#   Ω = Ω* + Ĉ T̂ Ĉ'  (Wall & Amemiya 2000, 2003)
# Otherwise, fall back to Ω = Ω* (Stage 2 only).
compute2SMMSE <- function(alpha_hat, M_hat_inv, R, f_eta, n,
                           C_hat = NULL, T_hat = NULL, ell = NULL) {
  # Residuals (always computed for variance estimation output)
  e_hat <- f_eta - R %*% alpha_hat

  # Omega* from ell contributions
  if (is.null(ell)) {
    # Fallback: raw ell = e_hat * R
    ell <- as.vector(e_hat) * R
  }
  Omega_star <- crossprod(ell) / n

  # Full sandwich meat
  if (!is.null(C_hat) && !is.null(T_hat)) {
    Omega <- Omega_star + C_hat %*% T_hat %*% t(C_hat)
  } else {
    Omega <- Omega_star
  }

  # V(alpha) = (1/n) M_hat^{-1} Omega M_hat^{-1}
  vcov_alpha <- (M_hat_inv %*% Omega %*% M_hat_inv) / n
  se <- sqrt(pmax(diag(vcov_alpha), 0))
  names(se) <- names(alpha_hat)

  list(vcov = vcov_alpha, se = se, Omega_star = Omega_star, Omega = Omega,
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


# ===========================================================================
# Full sandwich SE correction: Ω = Ω* + Ĉ T̂ Ĉ'
# (Wall & Amemiya 2000, 2003)
# ===========================================================================

# Get CFA parameters relevant to factor score computation (Lambda, nu, Theta).
# Psi (factor covariance) does not affect Bartlett factor scores, so we skip it.
# Returns a data.frame with columns: matrix, row, col, par_idx, est
get_cfa_relevant_params <- function(cfa_fit) {
  free_mats <- lavaan::lavInspect(cfa_fit, "free")
  est_mats  <- lavaan::lavInspect(cfa_fit, "est")

  # Matrices that affect factor scores: lambda, nu (intercepts), theta
  relevant_matrices <- c("lambda", "nu", "theta")

  params <- data.frame(
    matrix  = character(0),
    row     = integer(0),
    col     = integer(0),
    par_idx = integer(0),
    est     = numeric(0),
    stringsAsFactors = FALSE
  )

  for (mat_name in relevant_matrices) {
    if (!mat_name %in% names(free_mats)) next
    free_mat <- free_mats[[mat_name]]
    est_mat  <- est_mats[[mat_name]]

    # Find all free parameters (index > 0)
    idx <- which(free_mat > 0, arr.ind = TRUE)
    if (nrow(idx) == 0) next

    for (k in seq_len(nrow(idx))) {
      i <- idx[k, 1]
      j <- idx[k, 2]
      params <- rbind(params, data.frame(
        matrix  = mat_name,
        row     = i,
        col     = j,
        par_idx = free_mat[i, j],
        est     = est_mat[i, j],
        stringsAsFactors = FALSE
      ))
    }
  }

  # Sort by par_idx for consistent ordering
  params <- params[order(params$par_idx), ]
  rownames(params) <- NULL
  params
}


# Reconstruct factor scores from perturbed CFA parameters.
# param_info: data.frame from get_cfa_relevant_params()
# param_values: numeric vector of perturbed parameter values (same length as nrow(param_info))
# Z: n x p raw observed data matrix
# lvs: character vector of latent variable names (in desired order)
# base_Lambda, base_tau, base_Theta: baseline matrices (will be modified by perturbation)
reconstruct_stage1_from_perturbed <- function(param_info, param_values,
                                               Z, lvs,
                                               base_Lambda, base_tau, base_Theta) {
  Lambda <- base_Lambda
  tau    <- base_tau
  Theta  <- base_Theta

  # Apply perturbed values

  for (k in seq_len(nrow(param_info))) {
    mat <- param_info$matrix[k]
    i   <- param_info$row[k]
    j   <- param_info$col[k]
    val <- param_values[k]

    if (mat == "lambda") {
      Lambda[i, j] <- val
    } else if (mat == "nu") {
      tau[i] <- val
    } else if (mat == "theta") {
      Theta[i, j] <- val
      Theta[j, i] <- val  # symmetric
    }
  }

  # Recompute Bartlett factor scores
  Theta_inv <- solve(Theta)
  LtTiL     <- crossprod(Lambda, Theta_inv %*% Lambda)
  Sigma_ee  <- solve(LtTiL)
  W         <- Sigma_ee %*% crossprod(Lambda, Theta_inv)

  Z_centered <- sweep(Z, 2, tau)
  f_hat      <- Z_centered %*% t(W)
  colnames(f_hat) <- colnames(Lambda)

  # Reorder to match lvs
  f_hat    <- f_hat[, lvs, drop = FALSE]
  Sigma_ee <- Sigma_ee[lvs, lvs, drop = FALSE]
  W        <- W[lvs, , drop = FALSE]

  list(f_hat = f_hat, W = W, Sigma_ee = Sigma_ee, Theta = Theta)
}


# Compute observation-level estimating equation contributions.
# If error_moment_fn is provided, uses the exact bias-corrected formulation
# (Paper Eq 13-14): ℓ_t = m̂_t - M̂_t α̂
# Otherwise falls back to the raw approximation: ℓ_t = ê_t * r_t
# Returns an n x p_reg matrix of contributions.
compute_ell_contributions <- function(f_hat, alpha_hat, reg_multiidx,
                                       eta_name, factor_names,
                                       error_moment_fn = NULL) {
  if (!is.null(error_moment_fn)) {
    return(compute_bc_ell(f_hat, alpha_hat, reg_multiidx,
                          eta_name, factor_names, error_moment_fn))
  }

  n     <- nrow(f_hat)
  n_reg <- length(reg_multiidx)

  # Build regressor matrix R
  R <- matrix(0, nrow = n, ncol = n_reg)
  R[, 1] <- 1  # intercept
  for (j in seq_len(n_reg)[-1]) {
    r_j    <- reg_multiidx[[j]]
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

  # Residuals and contributions
  f_eta <- f_hat[, eta_name]
  e_hat <- f_eta - R %*% alpha_hat
  ell   <- as.vector(e_hat) * R  # n x n_reg

  ell
}


# Compute Ĉ matrix via central finite differences.
# Ĉ[j, m] = (1/n) Σ_t [ℓ_t^j(θ+h·e_m) - ℓ_t^j(θ-h·e_m)] / (2h)
#
# Returns a p_reg x M matrix where p_reg = number of regressors,
# M = number of relevant CFA parameters.
compute_C_hat <- function(stage1, alpha_hat, reg_multiidx,
                           eta_name, factor_names, lvs,
                           param_info, h_scale = 1e-5,
                           eps_cumulants = NULL) {
  n     <- stage1$n
  M     <- nrow(param_info)
  n_reg <- length(reg_multiidx)

  base_values <- param_info$est
  C_hat <- matrix(0, nrow = n_reg, ncol = M)

  for (m in seq_len(M)) {
    h <- max(abs(base_values[m]), 1) * h_scale

    # Perturb +h
    vals_plus     <- base_values
    vals_plus[m]  <- vals_plus[m] + h

    stage1_plus <- reconstruct_stage1_from_perturbed(
      param_info, vals_plus, stage1$Z, lvs,
      stage1$Lambda, stage1$tau, stage1$Theta
    )

    # Build perturbed error_moment_fn if eps_cumulants provided
    emf_plus <- NULL
    if (!is.null(eps_cumulants)) {
      eps_cum_plus <- eps_cumulants
      eps_cum_plus[[2]] <- diag(stage1_plus$Theta)
      emf_plus <- make_error_moment_fn(stage1_plus$W, eps_cum_plus)
    }

    ell_plus <- compute_ell_contributions(
      stage1_plus$f_hat, alpha_hat, reg_multiidx, eta_name, factor_names,
      error_moment_fn = emf_plus
    )

    # Perturb -h
    vals_minus     <- base_values
    vals_minus[m]  <- vals_minus[m] - h

    stage1_minus <- reconstruct_stage1_from_perturbed(
      param_info, vals_minus, stage1$Z, lvs,
      stage1$Lambda, stage1$tau, stage1$Theta
    )

    emf_minus <- NULL
    if (!is.null(eps_cumulants)) {
      eps_cum_minus <- eps_cumulants
      eps_cum_minus[[2]] <- diag(stage1_minus$Theta)
      emf_minus <- make_error_moment_fn(stage1_minus$W, eps_cum_minus)
    }

    ell_minus <- compute_ell_contributions(
      stage1_minus$f_hat, alpha_hat, reg_multiidx, eta_name, factor_names,
      error_moment_fn = emf_minus
    )

    # Central difference: Ĉ_m = (1/n) Σ_t [ℓ_t(θ+h) - ℓ_t(θ-h)] / (2h)
    C_hat[, m] <- colMeans(ell_plus - ell_minus) / (2 * h)
  }

  C_hat
}

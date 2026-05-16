# to do: add comments explaining everything and probably connecting to papers (2000; 2003)
# further check the analytical SEs for the unspecified distribution case for var_e(e)
# prpbably also compare to LSAM as they should be relatively identical except when e's are non-normal

library(lavaan); library(MASS)

quad2smm_true_params <- function() {
  reliab <- 0.75
  ph1    <- 1
  mf1    <- 3
  gam0   <- 3
  gam1   <- 10
  gam2   <- -2
  b10    <- 2
  b20    <- 3
  b21    <- 0.8
  b31    <- 0.6
  phf2   <- 11 # hardcoded in the SAS file

  list(
    reliab = reliab, ph1 = ph1, mf1 = mf1,
    gam0 = gam0, gam1 = gam1, gam2 = gam2,
    b10 = b10, b20 = b20, b21 = b21, b31 = b31,
    ps1 = sqrt(ph1 * (1 - reliab) / reliab),
    ps2 = sqrt(b21^2 * ph1 * (1 - reliab) / reliab),
    ps3 = sqrt(b31^2 * ph1 * (1 - reliab) / reliab),
    ps4 = sqrt(phf2 * (1 - reliab) / reliab)
  )
}

quad2smm_simulate <- function(n, p) {
  y1u <- runif(n)
  y1  <- sqrt(12) * (y1u - 0.5) # uniform on [-sqrt(3), sqrt(3)]
  x   <- matrix(rnorm(4 * n), n, 4)

  f1 <- p$mf1 + p$ph1 * y1
  Y1 <- f1 + p$ps1 * x[, 2]
  Y2 <- p$b10 + p$b21 * f1 + p$ps2 * x[, 3]
  Y3 <- p$b20 + p$b31 * f1 + p$ps3 * x[, 4]
  Y4 <- p$gam0 + p$gam1 * f1 + p$gam2 * f1^2 + p$ps4 * x[, 1]

  data.frame(y1 = Y1, y2 = Y2, y3 = Y3, y4 = Y4)
}

quad2smm_fit_meas <- function(Data) {
  model <- '
    f1 =~ 1*y1 + b21*y2 + b31*y3 + b41*y4
    y1 ~ m1*1
    y2 ~ b10*1
    y3 ~ b20*1
    y4 ~ b30*1
    y1 ~~ th1*y1
    y2 ~~ th2*y2
    y3 ~~ th3*y3
    y4 ~~ th4*y4
    f1 ~~ ph1*f1
    f1 ~ 0*1
    m2 := b10 - b21 * m1
    m3 := b20 - b31 * m1
  '
  lavaan::sem(model, data = Data, meanstructure = TRUE,
              estimator = "ML", baseline = FALSE)
}

quad2smm_extract_meas <- function(fit) {
  est <- lavaan::coef(fit)
  V   <- lavaan::vcov(fit)

  m1  <- unname(est["m1"])
  b21 <- unname(est["b21"])
  b31 <- unname(est["b31"])
  b10 <- unname(est["b10"])
  b20 <- unname(est["b20"])

  target <- c("m2", "m3", "b21", "b31", "th1", "th2", "th3")
  J <- matrix(0, length(target), length(est),
              dimnames = list(target, names(est)))
  J["m2", c("m1", "b10", "b21")] <- c(-b21, 1, -m1)
  J["m3", c("m1", "b20", "b31")] <- c(-b31, 1, -m1)
  J["b21", "b21"] <- 1
  J["b31", "b31"] <- 1
  J["th1", "th1"] <- 1
  J["th2", "th2"] <- 1
  J["th3", "th3"] <- 1
  cov7 <- J %*% V %*% t(J)

  list(
    b21  = b21,
    b31  = b31,
    th1  = unname(est["th1"]),
    th2  = unname(est["th2"]),
    th3  = unname(est["th3"]),
    ph1  = unname(est["ph1"]),
    m1   = m1,
    m2   = b10 - b21 * m1,
    m3   = b20 - b31 * m1,
    cov7 = cov7
  )
}

quad2smm_structural <- function(Data, mpar, n) {
  x    <- Data$y1
  z    <- Data$y2
  y    <- Data$y3
  newy <- Data$y4

  a1hat    <- mpar$b21
  c1hat    <- mpar$b31
  suuhat   <- mpar$th1
  saahat   <- mpar$th2
  seehat   <- mpar$th3
  covarian <- mpar$cov7

  mx <- mean(x); my <- mean(y); mz <- mean(z)
  mean1 <- mx
  mean2 <- mz
  mean3 <- my

  a0hat <- mz - a1hat * mx
  c0hat <- my - c1hat * mx

  out_names <- c("a1hat", "c1hat", "psiuu", "psiaa", "psiee", "seehat_dup",
                 "gamn0", "gamn1", "gamn2",
                 "segamn0", "segamn1", "segamn2", "sr2hat",
                 "bhatn80", "bhatn81", "bhatn82",
                 "bhat0", "bhat1", "bhat2",
                 "segam0_rob", "segam1_rob", "segam2_rob",
                 "bhat80", "bhat81", "bhat82",
                 "mean1", "mean2", "mean3")

  if (any(c(suuhat, saahat, seehat) < 5e-4)) {
    fhat <- if (suuhat < 5e-4)        x
            else if (saahat < 5e-4)   (z - a0hat) / a1hat
            else                      (y - c0hat) / c1hat
    hh   <- cbind(1, fhat, fhat^2)
    bhat <- drop(MASS::ginv(crossprod(hh)) %*% crossprod(hh, newy))
    out <- c(a1hat, c1hat, suuhat, saahat, seehat, seehat,
             bhat[1], bhat[2], bhat[3],
             NA, NA, NA, 0,
             bhat[1], bhat[2], bhat[3],
             bhat[1], bhat[2], bhat[3],
             NA, NA, NA,
             bhat[1], bhat[2], bhat[3],
             mean1, mean2, mean3)
    names(out) <- out_names
    return(out)
  }

  cc <- 1 / (suuhat * seehat * a1hat^2 +
             suuhat * saahat * c1hat^2 +
             saahat * seehat)

  sr2hat <- cc^2 * ((a1hat * seehat * suuhat)^2 * saahat +
                    (saahat * suuhat * c1hat)^2 * seehat +
                    (saahat * seehat)^2 * suuhat)

  fhat <- cc * (seehat * suuhat * a1hat * (z - a0hat) +
                saahat * suuhat * c1hat * (y - c0hat) +
                saahat * seehat * x)
  fhatsq <- fhat^2

  hh     <- cbind(1, fhat, fhatsq - sr2hat)
  mff    <- mean(fhatsq)
  mff2   <- mean(fhatsq - sr2hat)
  Mhat   <- crossprod(hh) / n
  mhyhat <- as.numeric(crossprod(hh, newy)) / n
  meanf  <- mean(fhat)

  Sigman <- matrix(0, 3, 3)
  Sigman[2, 2] <- sr2hat
  Sigman[2, 3] <- 2 * meanf * sr2hat
  Sigman[3, 2] <- Sigman[2, 3]
  Sigman[3, 3] <- 4 * sr2hat * mff - 2 * sr2hat^2

  bhatn <- drop(MASS::ginv(Mhat - Sigman) %*% mhyhat)

  v1 <- z - a0hat - a1hat * x
  v2 <- y - c0hat - c1hat * x

  dep3 <- c(mean(v1^3), mean(v1^2 * v2),
            mean(v1 * v2^2), mean(v2^3))
  projx3 <- matrix(0, 4, 3)
  projx3[1, 1] <-  1
  projx3[1, 2] <- -a1hat^3
  projx3[2, 2] <- -a1hat^2 * c1hat
  projx3[3, 2] <- -c1hat^2 * a1hat
  projx3[4, 2] <- -c1hat^3
  projx3[4, 3] <-  1
  ex3 <- drop(MASS::ginv(crossprod(projx3)) %*% crossprod(projx3, dep3))
  sa3hat <- ex3[1]; su3hat <- ex3[2]; se3hat <- ex3[3]

  dep4 <- c(mean(v1^4), mean(v1^3 * v2), mean(v1^2 * v2^2),
            mean(v2^3 * v1), mean(v2^4))
  inter4 <- c(
    6 * a1hat^2 * saahat * suuhat,
    3 * a1hat * c1hat * saahat * suuhat,
    saahat * seehat + a1hat^2 * suuhat * seehat + c1hat^2 * saahat * suuhat,
    3 * a1hat * c1hat * seehat * suuhat,
    6 * c1hat^2 * seehat * suuhat
  )
  projx4 <- matrix(0, 5, 3)
  projx4[1, 1] <- 1
  projx4[1, 2] <- a1hat^4
  projx4[2, 2] <- a1hat^3 * c1hat
  projx4[3, 2] <- a1hat^2 * c1hat^2
  projx4[4, 2] <- c1hat^3 * a1hat
  projx4[5, 2] <- c1hat^4
  projx4[5, 3] <- 1
  ex4 <- drop(MASS::ginv(crossprod(projx4)) %*%
                crossprod(projx4, dep4 - inter4))
  sa4hat <- max(sa3hat^2 / saahat, saahat^2, ex4[1])
  su4hat <- max(su3hat^2 / suuhat, suuhat^2, ex4[2])
  se4hat <- max(se3hat^2 / seehat, seehat^2, ex4[3])

  sr3hat <- cc^3 * ((a1hat * seehat * suuhat)^3 * sa3hat +
                    (saahat * suuhat * c1hat)^3 * se3hat +
                    (saahat * seehat)^3 * su3hat)

  sr4hat <- cc^4 * (
    (a1hat * seehat * suuhat)^4 * sa4hat +
    (saahat * suuhat * c1hat)^4 * se4hat +
    (saahat * seehat)^4 * su4hat +
    6 * seehat^2 * saahat^2 * suuhat^4 * a1hat^2 * c1hat^2 * saahat * seehat +
    6 * saahat^2 * suuhat^2 * seehat^4 * a1hat^2 * saahat * suuhat +
    6 * suuhat^2 * seehat^2 * saahat^4 * c1hat^2 * seehat * suuhat
  )

  Sigma <- matrix(0, 3, 3)
  Sigma[2, 2] <- sr2hat
  Sigma[2, 3] <- 2 * meanf * sr2hat + sr3hat
  Sigma[3, 2] <- Sigma[2, 3]
  Sigma[3, 3] <- 4 * mff2 * sr2hat + sr4hat + 4 * meanf * sr3hat - sr2hat^2

  bhat <- drop(MASS::ginv(Mhat - Sigma) %*% mhyhat)

  littlem <- t(hh * newy)
  subn    <- matrix(0, 3, n)
  for (t in seq_len(n)) {
    sigmat <- matrix(0, 3, 3)
    sigmat[2, 2] <- sr2hat
    sigmat[2, 3] <- 2 * fhat[t] * sr2hat
    sigmat[3, 2] <- sigmat[2, 3]
    sigmat[3, 3] <- 4 * sr2hat * fhatsq[t] - 2 * sr2hat^2
    bighn <- tcrossprod(hh[t, ]) - sigmat
    subn[, t] <- bighn %*% bhatn
  }
  dtn <- littlem - subn

  part1 <- numeric(3); part1[3] <- -mean(newy)
  part2 <- matrix(0, 3, 3)
  part2[1, 3] <- -1
  part2[2, 2] <- -1; part2[2, 3] <- -3 * meanf
  part2[3, 1] <- -1; part2[3, 2] <- -3 * meanf
  part2[3, 3] <- -6 * mff + 6 * sr2hat
  dldsig <- as.numeric(part1 - part2 %*% bhatn)

  sigvv <- matrix(c(
    saahat + a1hat^2 * suuhat, a1hat * c1hat * suuhat,
    a1hat * c1hat * suuhat,    seehat + c1hat^2 * suuhat
  ), 2, 2)
  invsig   <- solve(sigvv)
  modbet   <- matrix(c(a1hat, c1hat), nrow = 1)
  outer_bb <- crossprod(modbet)

  dsigth1 <- as.numeric(
    1
    - 2 * suuhat * drop(modbet %*% invsig %*% t(modbet))
    + suuhat^2   * drop(modbet %*% invsig %*% outer_bb %*% invsig %*% t(modbet))
  )
  dsigth2 <- suuhat^2 * drop(modbet %*% invsig %*% diag(c(1, 0)) %*%
                              t(invsig) %*% t(modbet))
  dsigth3 <- suuhat^2 * drop(modbet %*% invsig %*% diag(c(0, 1)) %*%
                              t(invsig) %*% t(modbet))

  e1  <- matrix(c(1, 0), 1, 2); e1c <- t(e1)
  e2v <- matrix(c(0, 1), 1, 2); e2c <- t(e2v)

  dsigb11 <- as.numeric(
    - suuhat^2 * drop(e1     %*% invsig %*% t(modbet))
    - suuhat^2 * drop(modbet %*% invsig %*% e1c)
    + suuhat^3 * drop(modbet %*% invsig %*%
                       (e1c %*% modbet + t(modbet) %*% e1) %*%
                       invsig %*% t(modbet))
  )
  dsigb12 <- as.numeric(
    - suuhat^2 * drop(e2v    %*% invsig %*% t(modbet))
    - suuhat^2 * drop(modbet %*% invsig %*% e2c)
    + suuhat^3 * drop(modbet %*% invsig %*%
                       (e2c %*% modbet + t(modbet) %*% e2v) %*%
                       invsig %*% t(modbet))
  )

  dsigdthe <- c(0, 0, dsigb11, dsigb12, dsigth1, dsigth2, dsigth3)
  coeff    <- rbind(cbind(diag(4), matrix(0, 4, 3)), dsigdthe)
  covnew   <- coeff %*% covarian %*% t(coeff)

  biggam <- -suuhat * matrix(c(a1hat, c1hat), nrow = 1) %*%
    solve(diag(c(saahat, seehat)) +
          suuhat * crossprod(matrix(c(a1hat, c1hat), nrow = 1)))

  average <- matrix(0, 3, 4)
  for (t in seq_len(n)) {
    p1 <- numeric(3); p1[2] <- newy[t]; p1[3] <- 2 * fhat[t] * newy[t]
    p2 <- matrix(0, 3, 3)
    p2[1, 2] <- 1
    p2[1, 3] <- 2 * fhat[t]
    p2[2, 1] <- 1
    p2[2, 2] <- 2 * fhat[t]
    p2[2, 3] <- 3 * fhatsq[t] - 3 * sr2hat
    p2[3, 1] <- 2 * fhat[t]
    p2[3, 2] <- 3 * fhatsq[t] - 3 * sr2hat
    p2[3, 3] <- 4 * fhatsq[t] * fhat[t] - 12 * fhat[t] * sr2hat
    dldf <- p1 - p2 %*% bhatn

    betat <- cbind(biggam, biggam * x[t])
    average <- average + (dldf %*% betat) / n
  }

  climit   <- cbind(average, dldsig)
  addition <- climit %*% covnew %*% t(climit)
  omegahn  <- tcrossprod(dtn) / n^2 + addition
  Minv     <- MASS::ginv(Mhat - Sigman)
  vhatn    <- Minv %*% omegahn %*% Minv

  segamn0 <- sqrt(max(0, vhatn[1, 1]))
  segamn1 <- sqrt(max(0, vhatn[2, 2]))
  segamn2 <- sqrt(max(0, vhatn[3, 3]))

  subn_rob <- matrix(0, 3, n)
  for (t in seq_len(n)) {
    sigmat <- matrix(0, 3, 3)
    sigmat[2, 2] <- sr2hat
    sigmat[2, 3] <- 2 * fhat[t] * sr2hat + sr3hat
    sigmat[3, 2] <- sigmat[2, 3]
    sigmat[3, 3] <- 4 * sr2hat * fhatsq[t] + sr4hat +
                    4 * sr3hat * fhat[t] - 5 * sr2hat^2
    bigh_rob <- tcrossprod(hh[t, ]) - sigmat
    subn_rob[, t] <- bigh_rob %*% bhat
  }
  dt_rob <- littlem - subn_rob

  part2_sig_rob <- matrix(0, 3, 3)
  part2_sig_rob[1, 3] <- -1
  part2_sig_rob[2, 2] <- -1
  part2_sig_rob[2, 3] <- -3 * meanf
  part2_sig_rob[3, 1] <- -1
  part2_sig_rob[3, 2] <- -3 * meanf
  part2_sig_rob[3, 3] <- -6 * mff + 12 * sr2hat
  dldsig_rob <- as.numeric(part1 - part2_sig_rob %*% bhat)

  average_rob <- matrix(0, 3, 4)
  for (t in seq_len(n)) {
    p1 <- numeric(3); p1[2] <- newy[t]; p1[3] <- 2 * fhat[t] * newy[t]
    p2 <- matrix(0, 3, 3)
    p2[1, 2] <- 1
    p2[1, 3] <- 2 * fhat[t]
    p2[2, 1] <- 1
    p2[2, 2] <- 2 * fhat[t]
    p2[2, 3] <- 3 * fhatsq[t] - 3 * sr2hat
    p2[3, 1] <- 2 * fhat[t]
    p2[3, 2] <- 3 * fhatsq[t] - 3 * sr2hat
    p2[3, 3] <- 4 * fhatsq[t] * fhat[t] - 12 * fhat[t] * sr2hat - 4 * sr3hat
    dldf <- p1 - p2 %*% bhat

    betat <- cbind(biggam, biggam * x[t])
    average_rob <- average_rob + (dldf %*% betat) / n
  }

  climit_rob   <- cbind(average_rob, dldsig_rob)
  addition_rob <- climit_rob %*% covnew %*% t(climit_rob)
  omega_rob    <- tcrossprod(dt_rob) / n^2 + addition_rob
  Minv_rob     <- MASS::ginv(Mhat - Sigma)
  vhat_rob     <- Minv_rob %*% omega_rob %*% Minv_rob

  segam0_rob <- sqrt(max(0, vhat_rob[1, 1]))
  segam1_rob <- sqrt(max(0, vhat_rob[2, 2]))
  segam2_rob <- sqrt(max(0, vhat_rob[3, 3]))

  myy <- sum(newy^2) / n
  cm  <- rbind(c(myy, mhyhat), cbind(mhyhat, Mhat))
  bias_adjusted <- function(Sig) {
    sigma0 <- rbind(rep(0, 4), cbind(rep(0, 3), Sig))
    rcm <- tryCatch(t(solve(chol(cm))), error = function(e) NULL)
    if (is.null(rcm)) return(rep(NA_real_, 3))
    msm <- rcm %*% sigma0 %*% t(rcm)
    ll  <- Re(eigen(msm, only.values = TRUE)$values)
    l   <- 1 / max(ll)
    lk  <- min(1, l - 1 / n)
    a8  <- lk - 8 / n
    drop(MASS::ginv(Mhat - Sig * a8) %*% mhyhat)
  }
  bhatn8 <- bias_adjusted(Sigman)
  bhat8  <- bias_adjusted(Sigma)

  out <- c(a1hat, c1hat, suuhat, saahat, seehat, seehat,
           bhatn[1], bhatn[2], bhatn[3],
           segamn0, segamn1, segamn2, sr2hat,
           bhatn8[1], bhatn8[2], bhatn8[3],
           bhat[1], bhat[2], bhat[3],
           segam0_rob, segam1_rob, segam2_rob,
           bhat8[1], bhat8[2], bhat8[3],
           mean1, mean2, mean3)
  names(out) <- out_names
  out
}

quad2smm_run <- function(n_reps = 1000, n = 500, seed = 1234,
                          verbose = TRUE) {
  set.seed(seed)
  p <- quad2smm_true_params()

  col_names <- c("a1hat", "c1hat", "psiuu", "psiaa", "psiee", "seehat_dup",
                 "gamn0", "gamn1", "gamn2",
                 "segamn0", "segamn1", "segamn2", "sr2hat",
                 "bhatn80", "bhatn81", "bhatn82",
                 "bhat0", "bhat1", "bhat2",
                 "segam0_rob", "segam1_rob", "segam2_rob",
                 "bhat80", "bhat81", "bhat82",
                 "mean1", "mean2", "mean3")
  out <- matrix(NA_real_, n_reps, length(col_names),
                dimnames = list(NULL, col_names))

  for (r in seq_len(n_reps)) {
    dat <- quad2smm_simulate(n, p)
    fit <- tryCatch(quad2smm_fit_meas(dat), error = function(e) NULL)
    if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) {
      if (verbose) cat(sprintf("rep %d: measurement model failed\n", r))
      next
    }
    mpar <- quad2smm_extract_meas(fit)
    res  <- tryCatch(quad2smm_structural(dat, mpar, n),
                     error = function(e) {
                       if (verbose) cat(sprintf("rep %d: structural failed: %s\n",
                                                r, conditionMessage(e)))
                       rep(NA_real_, length(col_names))
                     })
    out[r, ] <- res
    if (verbose && r %% 50 == 0)
      cat(sprintf("rep %4d / %d\n", r, n_reps))
  }
  as.data.frame(out)
}

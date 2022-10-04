#' @importFrom MGLM kr
#' @importFrom prodlim row.match
#'
fit.gplsvcm.worker <- function(iter.Tr, d, r, n.layer, V, Tr,
                               Y.all, X.all, Z.all, S.all, family, se = TRUE) {
  TV <- tdata(V, Tr)$TV
  dat.list <- subdata(iter.Tr, V, Tr, TV, n.layer, Y.all, X.all, Z.all, S.all)

  Vs <- dat.list$Vs
  Trs <- dat.list$Trs
  idxs <- dat.list$obs.ind
  t0 <- dat.list$t0
  Tr.set <- dat.list$subdomain.tri.set

  eta.cov <- NULL
  ZZ.tilde <- NULL
  eta <- NULL
  Var.fit <- NULL
  sigma <- NULL

  Ys <- Y.all[idxs]
  Xs <- X.all[idxs, ]
  Zs <- Z.all[idxs, ]
  Ss <- S.all[idxs, ]
  b.mtx <- basis(Vs, Trs, d = -1, r = 1, Ss)$Bi

  nq <- (d + 1) * (d + 2) / 2

  # Model fitting: Partially Linear Spatially Varying Coefficient Model
  if (!is.null(Xs) & !is.null(Zs)) {
    mfit <- mfit.gplsvcm(
      V = Vs, Tr = Trs, d, r, Y = Ys,
      X = Xs, Z = Zs, S = Ss, family
    )

    # extract basis coefficient w.r.t triangle of interest
    gammas <- mfit$gamma
    I <- matrix(1, nrow(Trs), ncol(Trs))
    tmp01 <- cbind(t0[1] * I[, 1], t0[2] * I[, 2], t0[3] * I[, 3])
    j <- which(apply(abs(Trs - tmp01), 1, sum) == 0)
    gamma.star <- gammas[((j - 1) * nq + 1):(j * nq), ]

    # linear coefficients
    eta <- mfit$eta
    eta.cov <- mfit$eta.cov

    # calculate elements in sd of eta
    if (se == TRUE) {
      Z.tilde <- mfit$Z.tilde
      Z.tilde <- kr(as.matrix(b.mtx), Z.tilde)
      # calculate V(mu.i) or g^{-1}'(\eta.i); two terms are equivalent
      Var.fit <- family$variance(mfit$fitted.values)
      # obtain sigma^2
      sigma <- mfit$mfit$sig2
      # ZZ.tilde
      ZZ.tilde <- crossprod(Z.tilde / sigma, Z.tilde)
    }
  }

  # spatially varying coefficient model ----
  if (!is.null(Xs) & is.null(Zs)) {
    mfit <- mfit.gsvcm(
      V = Vs, Tr = Trs, d, r, Y = Ys,
      X = Xs, S = Ss, family
    )
    gammas <- mfit$gamma
    I <- matrix(1, nrow(Trs), ncol(Trs))
    tmp01 <- cbind(t0[1] * I[, 1], t0[2] * I[, 2], t0[3] * I[, 3])
    j <- which(apply(abs(Trs - tmp01), 1, sum) == 0)
    gamma.star <- gammas[((j - 1) * nq + 1):(j * nq), ]
  }

  # results ----
  list(
    gamma.star = gamma.star,
    eta = eta,
    eta.cov = eta.cov,
    ZZ.tilde = ZZ.tilde,
    Var.fit = Var.fit,
    sigma = sigma,
    Tr.set = Tr.set
  )
}

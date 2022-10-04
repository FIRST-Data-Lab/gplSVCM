#' gplSVCM: Generalized Partially Linear Spatially Varying Coefficients Models
#'
#' mfit.gplsvcm(): This function fits a GPLSVCM on a subregion.
#'
#' @import BPST
#' @import Matrix
#' @import pracma
#' @import MGLM
#'
#' @param V	The nV by two matrix of vertices of a triangulation, where nV is the
#' number of vertices. Each row is the coordinates for a vertex.
#' @param Tr The triangulation matrix of dimention nT by three, where nT is the
#'  number of triangles in the triangulation. Each row is the indices of vertices in V.
#' @param d The degree of piecewise polynomials – default is 5, and usually d
#' is greater than one. -1 represents piecewise constant.
#' @param r	The smoothness parameter – default is 1, and 0 ≤ r < d.
#' @param Y A vector of response variable.
#' @param S A matrix containing the locations of new points.
#' @param Z A matrix containing the values of model covariates with
#' constant effects.
#' @param X A matrix containing the values of model covariates with
#' spatially varying effects.
#' @param family This is a family object specifying the distribution and link
#' to use.
#'

mfit.gplsvcm <- function(V, Tr, d, r, Y, X, Z, S, family) {
  Ball <- basis(V, Tr, d, r, S)
  B <- Ball$B
  Q2 <- Ball$Q2
  P <- Ball$K
  BQ2 <- as.matrix(B %*% Q2)
  PQ2 <- as.matrix(Matrix::crossprod(Q2, P) %*% Q2)
  X <- data.matrix(X)
  n <- length(Y)
  p1 <- ncol(Z)
  p2 <- ncol(X)
  nq <- ncol(Q2)

  dat.fit <- as.data.frame(cbind(Y, Z))
  pen.list <- list()

  for (ii in 1:p2) {
    dat.fit[[p1 + ii + 1]] <- kr(matrix(X[, ii], ncol = 1), BQ2)
    pen.list[[ii]] <- list(PQ2)
  }
  names(dat.fit)[1] <- "Y"
  names(dat.fit)[1 + 1:p1] <- paste0("Z", 1:p1)
  names(dat.fit)[1 + p1 + 1:p2] <- paste0("XB", 1:p2)
  names(pen.list) <- paste0("XB", 1:p2)
  formula <- as.formula(paste0("Y ~ 0 + ", paste0(names(dat.fit)[-1], collapse = " + ")))
  mfit <- gam(formula, family = family, paraPen = pen.list, data = dat.fit)

  eta <- mfit$coefficients[1:p1]
  theta <- matrix(mfit$coefficients[-(1:p1)], ncol = p2)
  gamma <- Q2 %*% theta
  beta <- BQ2 %*% theta
  eta.cov <- vcov(mfit, unconditional = TRUE)[1:p1, 1:p1]

  # calculate Z tilde
  Z.tilde <- NULL
  for (ii in 1:p2) {
    pen.list[[ii]] <- mfit$sp[ii] * PQ2
  }
  XB <- do.call("cbind", dat.fit[(2 + p1):(1 + p1 + p2)])
  fitted.values <- predict(mfit, dat.fit, type = "response")
  rho2 <- family$variance(fitted.values)
  XB <- as.vector(sqrt(rho2)) * XB
  Z <- as.vector(sqrt(rho2)) * Z
  P <- bdiag(pen.list)
  Z <- as.matrix(Z)
  Z.tilde <- Z - XB %*% solve(crossprod(XB), crossprod(XB, Z))

  list(
    mfit = mfit, eta = eta, eta.cov = eta.cov,
    theta = theta, gamma = gamma, beta = beta,
    Z.tilde = Z.tilde, fitted.values = fitted.values
  )
}

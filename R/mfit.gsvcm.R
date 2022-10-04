#' gplSVCM: Generalized Partially Linear Spatially Varying Coefficients Models
#'
#' mfit.gsvcm(): This function fits a generalized SVCM (GSVCM) on a subregion.
#'
#' @import BPST
#' @import Matrix
#' @import pracma
#' @import MGLM
#' @param V	The nV by two matrix of vertices of a triangulation, where nV is the
#' number of vertices. Each row is the coordinates for a vertex.
#' @param Tr The triangulation matrix of dimention nT by three, where nT is the
#'  number of triangles in the triangulation. Each row is the indices of vertices in V.
#' @param d The degree of piecewise polynomials – default is 5, and usually d
#' is greater than one. -1 represents piecewise constant.
#' @param r	The smoothness parameter – default is 1, and 0 ≤ r < d.
#' @param Y A vector of response variable.
#' @param S A matrix containing the locations of new points.
#' @param X A matrix containing the values of model covariates with
#' spatially varying effects.
#' @param family This is a family object specifying the distribution and link
#' to use.
#'

mfit.gsvcm <- function(V, Tr, d, r, Y, X, S, family){
  Ball = basis(V, Tr, d, r, S)
  B = Ball$B; dim(B)
  Q2 = Ball$Q2; dim(Q2)
  P = Ball$K; dim(P)
  BQ2 = as.matrix(B %*% Q2); dim(BQ2)
  PQ2 = as.matrix(Matrix::crossprod(Q2, P) %*% Q2); dim(PQ2)
  n = length(Y)
  p2 = ncol(X)
  nq = ncol(Q2)

  dat.fit = as.data.frame(Y)
  pen.list = list()
  for(ii in 1:p2){
    dat.fit[[ii + 1]] = kr(matrix(X[, ii], ncol = 1), BQ2)
    pen.list[[ii]] = list(PQ2)
  }
  names(dat.fit)[1] = "Y"
  names(dat.fit)[1 + 1:p2] = paste0("XB", 1:p2)
  names(pen.list) = paste0("XB", 1:p2)
  formula = as.formula(paste0("Y ~ 0 + ",
                              paste0(names(dat.fit)[-1], collapse = " + ")))
  mfit = gam(formula, family = family, paraPen = pen.list, data = dat.fit)

  theta = matrix(mfit$coefficients, ncol = p2)
  gamma = Q2 %*% theta
  beta = BQ2 %*% theta

  list(mfit = mfit, theta = theta, gamma = gamma,
       beta = beta, dat.fit = dat.fit)
}

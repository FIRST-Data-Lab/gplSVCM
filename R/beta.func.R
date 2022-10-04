#' gplSVCM: generalized partially linear model with spatially varying coefficients
#'
#' This function generates the values of true coefficient functions at
#' given points.
#'
#' @import BPST
#' @param S A matrix of the locations of points to be evaluated.
#' @param V	The nV by two matrix of vertices of a triangulation, where nV is the
#' number of vertices. Each row is the coordinates for a vertex.
#' @param Tr The triangulation matrix of dimention nT by three, where nT is the
#'  number of triangles in the triangulation. Each row is the indices of vertices in V.
#' @return
#' \item{beta0}{A matrix of the values of beta0 and beta1 at location points S.}
#' @export
beta.func <- function(S, V, Tr) {
  xx <- S[, 1]
  yy <- S[, 2]
  x0 <- unique(xx)
  y0 <- unique(yy)
  n <- nrow(S)
  ind.in <- which(insideVT(V, Tr, xx, yy))
  ind.all <- 1:n
  ind.out <- setdiff(ind.all, ind.in)

  r0 <- 0.1
  r <- 0.5
  l <- 3
  b <- 1
  q <- pi * r / 2
  a <- rep(0, n)
  d <- rep(0, n)

  # part 1
  ind <- (xx >= 0) & (yy > 0)
  a[ind] <- q + xx[ind]
  d[ind] <- yy[ind] - r
  # part 2
  ind <- (xx >= 0) & (yy <= 0)
  a[ind] <- -q - xx[ind]
  d[ind] <- -r - yy[ind]
  # part 3
  ind <- (xx < 0)
  a[ind] <- -atan(yy[ind] / xx[ind]) * r
  d[ind] <- sqrt(xx[ind]^2 + yy[ind]^2) - r

  f1 <- a * b + d^2

  ind <- (abs(d) > r - r0) | (xx > l & (xx - l)^2 + d^2 > (r - r0)^2)
  f1[ind] <- NA
  f1[ind.out] <- NA
  f2 <- -f1

  beta0 <- cbind(f1, f2)

  return(beta0)
}

#' gplSVCM: generalized partially linear model with spatially varying coefficients
#'
#' This function identifies points falling within a triangulation.
#'
#' @import BPST
#'
dat.dc <- function(Vs, Trs, Y, X, Z, S) {
  nT <- nrow(Trs)
  umin <- min(Vs[, 1])
  umax <- max(Vs[, 1])
  vmin <- min(Vs[, 2])
  vmax <- max(Vs[, 2])
  uu <- S[, 1]
  vv <- S[, 2]
  i <- which(uu >= umin & uu <= umax)
  j <- which(vv[i] >= vmin & vv[i] <= vmax)
  Sij <- cbind(uu[i[j]], vv[i[j]])
  uu <- Sij[, 1]
  vv <- Sij[, 2]
  Yij <- Y[i[j]]
  if (!is.null(X)) {
    if (ncol(X) > 1) {
      Xij <- X[i[j], ]
    }
  }
  Zij <- Z[i[j], ]
  Index0 <- i[j]

  Ys <- NULL
  Xs <- NULL
  Zs <- NULL
  Ss <- NULL
  tri.obs.size <- c()
  Index <- NULL
  for (ii in 1:nT) {
    lam <- BPST::bary(Vs[Trs[ii, 1], ], Vs[Trs[ii, 2], ], Vs[Trs[ii, 3], ], uu, vv)
    lam1 <- lam$lam1
    lam2 <- lam$lam2
    lam3 <- lam$lam3
    I <- which(lam1 >= 0 & lam2 >= 0 & lam3 >= 0)
    Ys <- c(Ys, Yij[I])
    if (!is.null(X)) {
      if (ncol(X) > 1) {
        Xs <- rbind(Xs, Xij[I, ])
      }
    }
    Zs <- rbind(Zs, Zij[I, ])
    Ss <- rbind(Ss, cbind(uu[I], vv[I]))
    tri.obs.size <- c(tri.obs.size, length(I))
    Index <- c(Index, Index0[I])
  }
  if (!is.null(X)) {
    if (ncol(X) == 1) {
      Xs <- matrix(1, ncol = 1, nrow = dim(Zs)[1])
    }
  }
  list(
    Ys = Ys, Xs = Xs, Zs = Zs, Ss = Ss, tri.obs.size = tri.obs.size,
    Index = Index
  )
}

#' @import Triangulation MGLM
#' @importFrom prodlim row.match
#' @importFrom Matrix rankMatrix
#'
subdata <- function(ii, V, Tr, TV, n.layer, Y, X, Z, S) {
  cat("Working on triangle:", ii, "\n")
  # obtain sub-triangulation
  Tr0 <- Tr[ii, ]
  VT <- ring(Tr0, V, Tr, TV, n.layer)
  Vs <- VT$V1
  Trs <- VT$T1
  t0 <- VT$Tr0

  # obtain data within each sub-domain
  dats <- dat.dc(Vs, Trs, Y, X = X, Z, S)
  Xs <- dats$Xs
  Zs <- dats$Zs
  Ss <- dats$Ss

  # check whether expand the sub-domain.
  # condition 1 whether the Xs is column-full rank matrix
  # condition 2 whether the Zs is column-full rank matrix
  # condition 3 whether there are enough observation points for the sub-domain
  # condition 4 check number of triangles within each sub-domain
  if(!is.null(Xs) & !is.null(Zs)) {
    conditions <- c(
      rankMatrix(Zs) == ncol(Zs), rankMatrix(Xs) == ncol(Xs),
      nrow(Ss) / nrow(Trs), nrow(Trs)
    )
  }

  if(is.null(Xs) & !is.null(Zs)) {
    conditions <- c(
      rankMatrix(Zs) == ncol(Zs),
      nrow(Ss) / nrow(Trs), nrow(Trs)
    )
  }

  if(!is.null(Xs) & is.null(Zs)) {
    conditions <- c(
      rankMatrix(Xs) == ncol(Xs),
      nrow(Ss) / nrow(Trs), nrow(Trs)
    )
  }

  # extract information within each sub-triangulation
  # index of triangles within sub-domain iter
  subdomain.tri.set <- row.match(data.frame(VT$T.sub), data.frame(Tr))
  # number of observations within each triangle
  tri.obs.size <- dats$tri.obs.size
  # index of obs. w.r.t the whole sample
  obs.ind <- dats$Index

  list(
    Vs = Vs, Trs = Trs, t0 = t0, subdomain.tri.set = subdomain.tri.set,
    tri.obs.size = tri.obs.size, obs.ind = obs.ind,
    conditions = conditions
  )
}

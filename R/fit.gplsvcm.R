#' gplSVCM: Generalized Partially Linear Spatially Varying Coefficients Models
#'
#' fit.gplsvcm(): This function fits a GPLSVCM based on divide-and-conquer and domain decomposition methods.
#'
#' @import parallel BPST
#' @import Triangulation
#' @importFrom Matrix sparseMatrix Matrix
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
#' @param n.layer Number of layers of triangles used to define the neighborhood.
#' @param se An indicator of whether to return the estimated variance-covariance matrix of estimator of the linear components.
#' @param ns Number of parallel cores to use.
#'
#' @examples
#' # Triangulation
#' data("V.2"); data("Tr.2");
#' V = V.2; Tr = Tr.2;
#' Triangulation::TriPlot(V, Tr)
#' # Parameters
#' d <- 2; r <- 1;
#' family.name <- "gau"
#' if (family.name == "bin") {
#' family <- binomial()
#' } else if (family.name == "gau") {
#' family <- gaussian()
#' }
#' n.layer = 2; ns <- 2;
#' n <- 10000;
#' dat <- dat.generator(n, V.2, Tr.2, family.name)
#' Y <- dat$y
#' X <- cbind(dat$x1, dat$x2)
#' Z <- cbind(dat$z1, dat$z2, dat$z3)
#' S <- cbind(dat$s1, dat$s2)
#' mfit <- fit.gplsvcm(V, Tr, d, r,
#' Y = Y, X = X, Z = Z, S = S,
#' family = family, n.layer = n.layer, se = TRUE, ns = ns)
#' # Results: eta.hat and its variance-covariance matrix
#' mfit$eta.hat
#' mfit$eta.hat.cov
#'
#' # Generate population grid points
#' uu.grid <- seq(-0.9, 3.4, 0.02)
#' vv.grid <- seq(-0.9, 0.9, 0.02)
#' S.grid <- as.matrix(expand.grid(uu.grid, vv.grid))
#' # plots of estimated coefficient functions
#' plots.gplsvcm(mfit, S.grid)
#'
#' @export
#'

fit.gplsvcm <- function(V, Tr, d, r, Y, S, Z, X, family,
                        n.layer = 2, se = FALSE, ns = 2) {
  this.call <- match.call()

  TV <- tdata(V, Tr)$TV
  n.T <- nrow(Tr)

  # Model fitting:
  # Step 1: Divide (Domain Decomposition)
  results.all <- mclapply(1:nrow(Tr),
    FUN = fit.gplsvcm.worker, mc.cores = ns,
    d = d, r = r, n.layer = n.layer,
    Y.all = Y, Z.all = Z, X.all = X, S.all = S,
    family = family, se = se, V = V, Tr = Tr,
    mc.preschedule = FALSE
  )

  gamma.star <- list()
  eta <- list()
  eta.cov <- list()
  Tr.set <- list()
  ZZ.tilde <- list()

  # combine results from workers ----
  for (argument in names(results.all[[1]])) {
    assign(argument, lapply(results.all, "[[", argument))
  }

  # Step 2: Aggregation
  # Step 2.1 Re-order and stack the spatially varying coefficients, gamma
  # Step 2.2 Adjust to guarantee the global smoothness
  if (!is.null(X)) {
    gamma.star.all <- do.call("rbind", gamma.star)
    if (ncol(X) == 1) gamma.star.all <- c(gamma.star.all)

    H <- as.matrix(smoothness(V, Tr, d, r))
    a <- H %*% gamma.star.all
    HH <- crossprod(t(H))
    nH <- nrow(HH)
    b <- chol2inv(chol(HH + 1e-12 * diag(nH))) %*% a
    gamma <- gamma.star.all - t(H) %*% b
  }

  # linear coefficients eta
  eta.hat <- NULL
  eta.cov.comb <- NULL
  eta.hat.cov <- NULL

  # combine eta
  if (!is.null(Z)) {
    eta.all.2 <- eta
    eta.cov.inV.all <- matrix(0, ncol(Z), ncol(Z))
    for (iter in 1:n.T) {
      eta.cov.inv <- solve(eta.cov[[iter]])
      eta.all.2[[iter]] <- eta.cov.inv %*% eta[[iter]]
      eta.cov.inV.all <- eta.cov.inV.all + eta.cov.inv
    }
    eta.cov.comb <- solve(eta.cov.inV.all)
    eta.all.2 <- do.call("cbind", eta.all.2)
    eta.hat <- eta.cov.comb %*% rowSums(eta.all.2)

    ZZ.all <- matrix(0, ncol(Z), ncol(Z))

    if (se == TRUE) {
      # calculate Z.i'D.iZ.i for each observation point
      ZZ.all <- mclapply(1:n.T,
        FUN = ZZ.func, mc.cores = ns, Tr.set = Tr.set,
        ZZ.tilde = ZZ.tilde, p = ncol(Z)
      )
      ZZ <- Reduce(`+`, ZZ.all)
      eta.hat.cov <- eta.cov.comb %*% ZZ %*% eta.cov.comb
    }
  }

  res <- list(
    gamma.hat = gamma,
    gamma.star = gamma.star.all,
    eta.hat = eta.hat,
    eta.hat.cov = eta.hat.cov,
    V = V,
    Tr = Tr,
    d = d,
    r = r,
    family = family
  )
  res$call <- this.call;
  class(res) <- "gplsvcm"
  return(res)
}

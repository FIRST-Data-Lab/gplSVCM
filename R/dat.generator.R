#' gplSVCM: generalized partially linear model with spatially varying coefficients
#'
#' This function generates generate a dataset based on simulation studies.
#'
#' @import BPST MASS
#' @importFrom stats as.formula rbinom rnorm runif vcov binomial
#' @param n Number of observation points to be generated.
#' @param V	The nV by two matrix of vertices of a triangulation, where nV is the
#' number of vertices. Each row is the coordinates for a vertex.
#' @param Tr The triangulation matrix of dimention nT by three, where nT is the
#'  number of triangles in the triangulation. Each row is the indices of vertices in V.
#' @param family.name A family object specifying the distribution and link to use.
#' @param sig Standard deviation of the random error in gaussian distribution.
#' @param seed An integer of random seed.
#' @return
#' \item{dat}{A list of model components in a simlulated dataset.}
#' @export
#'
dat.generator <- function(n, V, Tr, family.name = "gau", seed = 2022) {
  set.seed(seed)
  # check if the points are within triangulation
  uu <- runif(2 * n, -0.9, 3.4)
  vv <- runif(2 * n, -0.9, 0.9)
  S.all <- cbind(uu, vv)
  beta.sample <- beta.func(S.all, V, Tr)
  sample.ind <- sample(which(!is.na(beta.sample[, 1])), n)
  f1 <- beta.sample[sample.ind, 1]
  f2 <- beta.sample[sample.ind, 2]
  S <- cbind(uu[sample.ind], vv[sample.ind])

  x1 <- rep(1, n)
  x.cov <- matrix(c(
    1, 0.5, 0.25, 0.125, 0.5, 1, 0.5, 0.25,
    0.25, 0.5, 1, 0.5, 0.125, 0.25, 0.5, 1
  ),
  nrow = 4, ncol = 4
  )
  mu <- rep(0, 4)
  x.all <- mvrnorm(n, mu = mu, Sigma = x.cov)
  x2 <- x.all[, 1]
  z1 <- x.all[, 2]
  z2 <- x.all[, 3]
  z3 <- x.all[, 4]

  eta0 <- c(1, -1, 1.5)
  eta <- f1 * x1 + f2 * x2 + eta0[1] * z1 + eta0[2] * z2 + eta0[3] * z3

  if (family.name == "bin") {
    ind <- which(!is.na(eta))
    family <- binomial()
    mu <- family$linkinv(eta)
    y <- mu
    for (ii in ind) {
      y[ii] <- rbinom(1, size = 1, prob = mu[ii])
    }
  } else {
    sig = 5
    mu <- eta
    e <- rnorm(n, 0, sig)
    y <- mu + e
  }

  dat <- as.data.frame(cbind(y, x1, x2, z1, z2, z3, S, mu, f1, f2))
  colnames(dat) <- c(
    "y", "x1", "x2", "z1", "z2", "z3",
    "s1", "s2", "mu", "f1", "f2"
  )
  return(dat)
}

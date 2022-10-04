rm(list = ls())

# library(devtools)
# install_github("funstatpackages/gplSVCM")
library(gplSVCM)

######################################################################
data("V.2"); data("Tr.2")
V = V.2; Tr = Tr.2
Triangulation::TriPlot(V, Tr)

d <- 2; r <- 1

# parameters
family.name <- "gau"
if (family.name == "bin") {
  family <- binomial()
} else if (family.name == "gau") {
  family <- gaussian()
}

n.layer = 2; ns <- 2

# generate random dataset ----
n <- 10000
dat <- dat.generator(n, V, Tr, family.name)

Y <- dat$y
X <- cbind(dat$x1, dat$x2)
Z <- cbind(dat$z1, dat$z2, dat$z3)
S <- cbind(dat$s1, dat$s2)

t0 <- proc.time()
mfit <- fit.gplsvcm(V, Tr, d, r,
                    Y = Y, X = X, Z = Z, S = S,
                    family = family,
                    n.layer = n.layer,
                    se = TRUE,
                    ns = ns)
t1 <- proc.time() - t0
# Estimated constant effects and its variance
mfit$eta.hat
mfit$eta.hat.cov

# Predict coefficient function at grid points
# generate population grid points
uu.grid <- seq(-0.9, 3.4, 0.02)
vv.grid <- seq(-0.9, 0.9, 0.02)
S.grid <- as.matrix(expand.grid(uu.grid, vv.grid))

# prediction
Z.new <- matrix(runif(3 * nrow(S.grid)), ncol = 3)
X.new <- cbind(1, runif(nrow(S.grid)))
pred <- predict.gplsvcm(mfit,
  S.new = S.grid,
  Z.new = Z.new, X.new = X.new,
  type = c("beta.terms", "response")
)
beta.pred <- pred$beta.terms

# calculate MISE
beta.true <- beta.func(S.grid, V, Tr)

# the MISE for BPST is
apply((beta.pred - beta.true)^2, 2, mean, na.rm = TRUE)

# plots of estimated coefficient functions
plots.gplsvcm(mfit, S.grid = S.grid)

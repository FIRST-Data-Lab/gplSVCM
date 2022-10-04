#' gplSVCM: generalized partially linear model with spatially varying coefficients
#'
#' This function generates raster plots for the estimated coefficient functions.
#'
#' @import ggplot2 colorRamps
#' @param mfit Fitted ``GPL-SVCM" object
#' @param S.grid A matrix containing the locations of points to be evaluated.
#' @return
#' \item{p.all}{A list of rasters plots of estimated coefficient functions.}
#'
#' @export
#'
plots.gplsvcm <- function(mfit, S.grid = NULL) {
  if (is.null(S.grid)) {
    z1.min <- min(mfit$V[, 1]) - 0.1
    z1.max <- max(mfit$V[, 1]) + 0.1
    z2.min <- min(mfit$V[, 2]) - 0.1
    z2.max <- max(mfit$V[, 2]) + 0.1
    n1 <- 50
    n2 <- 50
    u1 <- seq(z1.min, z1.max, length.out = n1)
    v1 <- seq(z2.min, z2.max, length.out = n2)
    S.grid <- as.matrix(expand.grid(u1, v1))
  }

  beta.mtx <- predict.gplsvcm(mfit, S.new = S.grid, type = "beta.terms")$beta.terms
  p.all <- list()
  for (iter in 1:ncol(beta.mtx)) {
    data.fit <- data.frame(
      s1 = S.grid[, 1], s2 = S.grid[, 2],
      fitted = beta.mtx[, iter]
    )

    p <- ggplot(data = data.fit, aes(x = s1, y = s2)) +
      geom_raster(aes(fill = fitted)) +
      coord_fixed(ratio = 2) +
      scale_fill_gradientn(colours = matlab.like(104), na.value = "transparent") +
      theme_bw() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(fill = paste0("beta", iter - 1))
    p.all[[iter]] <- p
  }

  return(p.all)
}

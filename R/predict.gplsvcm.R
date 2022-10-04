#' gplSVCM: generalized partially linear model with spatially varying coefficients
#'
#' This function generates predictions of model terms for a given new
#' dataset based on the estimated model from ``GPL-SVCM" object.
#'
#' @import BPST mgcv
#' @param mfit Fitted ``GPL-SVCM" object
#' @param S.new A matrix containing the locations of new points.
#' @param Z.new A matrix containing the values of model covariates with
#' constant effects.
#' @param X.new A matrix containing the values of model covariates with
#' spatially varying effects.
#' @param type A vector containing types of results to output. When type = "link"
#' (default), the linear predictor is returned. When type = "beta.terms",
#' the values of estimated beta functions at locations \code{S.new} are
#' returned. When type = "response", predictions on the mean of response
#' are returned.
#' @return
#' \item{result}{A list of predictions for different model terms.}
#' @export
#'
predict.gplsvcm <- function(mfit, S.new, Z.new = NULL, X.new = NULL,
                            type = c("beta.terms", "response")) {
  gamma <- mfit$gamma.hat

  basis.new <- basis(mfit$V, mfit$Tr, mfit$d, mfit$r, S.new)
  B.new <- basis.new$B
  ind.new <- basis.new$Ind.inside

  nv <- ncol(gamma)
  beta.mtx <- matrix(NA, nrow(S.new), nv)
  beta.mtx[ind.new, ] <- as.matrix(B.new %*% gamma)
  colnames(beta.mtx) <- paste0("beta", 1:nv)

  result <- list()
  if ("beta.terms" %in% type) {
    result$"beta.terms" <- beta.mtx
  }

  if ("link" %in% type) {
    if(is.null(mfit$eta.c)) {
      result$"link" <- rowSums(X.new * beta.mtx)
    } else {
      result$"link" <- Z.new %*% mfit$eta.c + rowSums(X.new * beta.mtx)
    }
  }

  if ("response" %in% type) {
    if(is.null(mfit$eta.c)) {
      result$"response" <-
        mfit$family$linkinv(rowSums(X.new * beta.mtx))
    } else {
      result$"response" <-
        mfit$family$linkinv(Z.new %*% mfit$eta.c + rowSums(X.new * beta.mtx))
    }
  }

  return(result)
}

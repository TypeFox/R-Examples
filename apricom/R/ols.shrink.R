#' Estimation of a Shrinkage Factor for Linear Regression
#'
#' Estimate a shrinkage factor for shrinkage-after-estimation techniques,
#' with application to linear regression models.
#'
#' This is an accessory function that works together with \code{\link{bootval}}, \code{\link{splitval}},
#' \code{\link{kcrossval}} and \code{\link{loocval}} to estimate a shrinkage factor. For further details,
#' see References. This function should not be used directly, and instead should
#' be called via one of the aforementioned shrinkage-after-estimation functions.
#'
#' @param b 1 x \code{m} matrix of regression coefficients, derived by resampling
#'                       or sample splitting
#' @param dat a \code{p} x \code{m} data matrix, where the final column is a continuous
#'            outcome variable. This dataset acts as a "test set" or "validation set".
#' @param sdm the shrinkage design matrix. This determines the regression coefficients
#'       that will be involved in the shrinkage process.
#' @return the function returns a shrinkage factor.
#'
#' @note Currently, this function can only derive a single shrinkage factor for a given
#' model, and is unable to estimate (weighted) predictor-specific shrinkage factors.
#'
#' @examples
#'## Shrinkage design matrix examples for a model with an
#'## intercept and 4 predictors:

#'## 1. Uniform shrinkage (default design within apricomp).
#'    sdm1 <- matrix(c(0, rep(1, 4)), nrow = 1)
#'    print(sdm1)

#'## 2. Non-uniform shrinkage; 1 shrinkage factor applied only to the
#'##    first two predictors
#'    sdm2 <- matrix(c(0, 1, 1, 0, 0), nrow = 1)
#'    print(sdm2)
#'
#' @references Harrell, F. E. \emph{"Regression modeling strategies: with applications
#'              to linear models, logistic regression, and survival analysis."} \emph{Springer}, (2001).
#' @references Steyerberg, E. W. \emph{"Clinical Prediction Models", Springer} (2009)

ols.shrink <- function(b, dat, sdm) {

  dat <- as.matrix(dat)
  datX <- dat[, 1:dim(dat)[2] - 1]
  sdm <- t(sdm)

# Adjustment for non-uniform shrinkage designs
  sdm.adj <- apply(1 - sdm, 1, min)

# Calculate linear predictor and combine with (adjusted) outcome variable
  LPY <- cbind(datX %*% diag(as.vector(b)) %*% sdm, dat[, dim(dat)[2]] -
                 datX %*% diag(as.vector(b)) %*% sdm.adj)

# Calculate shrinkage factor (regress outcome on linear predictor)
  s <- ols.rgr(LPY)

  return(s)
}

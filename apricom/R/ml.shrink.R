#' Estimation of a Shrinkage Factor for Logistic Regression
#'
#' Estimate a shrinkage factor for shrinkage-after-estimation techniques,
#' with application to logistic regression models.
#'
#' This function works together with \code{\link{bootval}}, \code{\link{splitval}},
#' \code{\link{kcrossval}} and \code{\link{loocval}} to estimate a shrinkage factor. For further details,
#' see References. This function should not be used directly, and instead should
#' be called via one of the aforementioned shrinkage-after-estimation functions.
#'
#' @param b 1 x \code{m} matrix of regression coefficients, derived by resampling or
#'           sample-splitting
#' @param dat a \code{p} x \code{m} data matrix, where the final column is a
#'        binary outcome variable. This dataset acts as a "test set" or "validation set".
#' @return the function returns a single shrinkage factor
#' @note Currently, this function can only derive a single shrinkage factor for a given
#' model, and is unable to estimate (weighted) predictor-specific shrinkage factors.
#' @references Harrell, F. E. \emph{"Regression modeling strategies: with applications
#'              to linear models, logistic regression, and survival analysis."} \emph{Springer}, (2001).
#' @references Steyerberg, E. W. \emph{"Clinical Prediction Models", Springer} (2009)

ml.shrink<- function(b, dat) {

  dat <- as.matrix(dat)
  nc <- dim(dat)[2]
  datX <- dat[, 1:(nc - 1)]
  Y <- dat[, nc]

# Calculate linear predictor and combine with (adjusted) outcome variable
LP <- datX %*% b
LPY <- cbind(LP, Y)

# Calculate shrinkage factor (regress outcome on linear predictor)
s <- ml.rgr(LPY)

return(s)
}

# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Penalty parameter for sparse LTS regression
#' 
#' Use bivariate winsorization to estimate the smallest value of the penalty 
#' parameter for sparse least trimmed squares regression that sets all 
#' coefficients to zero.
#' 
#' The estimation procedure is inspired by the calculation of the respective 
#' penalty parameter in the first step of the classical LARS algorithm.
#' First, two-dimensional data blocks consisting of the response with each 
#' predictor variable are cleaned via bivariate winsorization.  For each block, 
#' the following computations are then performed.  If an intercept is included 
#' in the model, the cleaned response is centered and the corresponding cleaned 
#' predictor is centered and scaled to have unit norm.  Otherwise the variables 
#' are not centered, but the predictor is scaled to have unit norm.  Finally, 
#' the dot product of the response and the corresponding predictor is 
#' computed.  The largest absolute value of those dot products, rescaled to fit 
#' the parametrization of the sparse LTS definition, yields the estimate of the 
#' smallest penalty parameter that sets all coefficients to zero.
#' 
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param normalize  a logical indicating whether the winsorized predictor 
#' variables should be normalized to have unit \eqn{L_{2}}{L2} norm (the 
#' default is \code{TRUE}).
#' @param intercept  a logical indicating whether a constant term should be 
#' included in the model (the default is \code{TRUE}).
#' @param const  numeric; tuning constant to be used in univariate 
#' winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in bivariate
#' winsorization (defaults to 0.95).
#' @param tol  a small positive numeric value used to determine singularity 
#' issues in the computation of correlation estimates for bivariate 
#' winsorization (see \code{\link{corHuber}}).
#' @param eps  a small positive numeric value used to determine whether the 
#' robust scale estimate of a variable is too small (an effective zero).
#' @param \dots additional arguments to be passed to 
#' \code{\link[=standardize]{robStandardize}}.
#' 
#' @return A robust estimate of the smallest value of the penalty parameter for 
#' sparse LTS regression that sets all coefficients to zero.
#' 
#' @author Andreas Alfons
#' 
#' @references
#' Alfons, A., Croux, C. and Gelper, S. (2013) Sparse least trimmed squares 
#' regression for analyzing high-dimensional large data sets. \emph{The Annals 
#' of Applied Statistics}, \bold{7}(1), 226--248.
#' 
#' Efron, B., Hastie, T., Johnstone, I. and Tibshirani, R. (2004) Least angle 
#' regression. \emph{The Annals of Statistics}, \bold{32}(2), 407--499.
#' 
#' Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust linear model 
#' selection based on least angle regression. \emph{Journal of the American 
#' Statistical Association}, \bold{102}(480), 1289--1299.
#' 
#' @seealso \code{\link{sparseLTS}}, \code{\link{winsorize}}
#' 
#' @example inst/doc/examples/example-lambda0.R
#' 
#' @keywords utilities robust
#' 
#' @export

lambda0 <- function(x, y, normalize = TRUE, intercept = TRUE, const = 2, 
                    prob = 0.95, tol = .Machine$double.eps^0.5, 
                    eps = .Machine$double.eps, ...) {
  # initializations
  n <- length(y)
  x <- as.matrix(x)
  if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
  normalize <- isTRUE(normalize)
  intercept <- isTRUE(intercept)
  # standardize data
  y <- robStandardize(y, eps=eps, ...)
  centerY <- attr(y, "center")
  scaleY <- attr(y, "scale")
  if(scaleY <= eps) stop("scale of response is too small")
  x <- robStandardize(x, eps=eps, ...)
  centerX <- attr(x, "center")
  scaleX <- attr(x, "scale")
  # drop variables with too small a scale
  keep <- which(scaleX > eps)
  if(length(keep) == 0) stop("scale of all predictors is too small")
  x <- x[, keep, drop=FALSE]
  centerX <- centerX[keep]
  scaleX <- scaleX[keep]
  # compute largest lambda
  corY <- sapply(seq_len(ncol(x)), function(j) {
    # winsorize bivariate data and transform back to original scale
    tmp <- winsorize(cbind(x[, j], y), const=const, 
                     prob=prob, standardized=TRUE, tol=tol)
    xj <- tmp[, 1] * scaleX[j] + centerX[j]
    y <- tmp[, 2] * scaleY + centerY
    # in case of intercept, sweep out means of winsorized data
    if(intercept) {
      xj <- xj - mean(xj)
      y <- y - mean(y)
    }
    # normalize cleaned x variable with respect to L2 norm
    if(normalize) xj <- xj / sqrt(sum(xj^2))
    # return 
    drop(t(y) %*% xj)
  })
  max(abs(corY)) * 2 / nrow(x)
}

# ------------------------------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
#
# functionality for correlations via winsorization is based on 
# code by Jafar A. Khan, Stefan Van Aelst and Ruben H. Zamar
# ------------------------------------------------------------

#' Robust correlation based on winsorization.
#' 
#' Compute a robust correlation estimate based on winsorization, i.e., by 
#' shrinking outlying observations to the border of the main part of the data.
#' 
#' The borders of the main part of the data are defined on the scale of the 
#' robustly standardized data.  In univariate winsorization, the borders for 
#' each variable are given by \eqn{+/-}\code{const}, thus a symmetric 
#' distribution is assumed.  In adjusted univariate winsorization, the borders 
#' for the two diagonally opposing quadrants containing the minority of the 
#' data are shrunken by a factor that depends on the ratio between the number of 
#' observations in the major and minor quadrants.  It is thus possible to 
#' better account for the bivariate structure of the data while maintaining 
#' fast computation.  In bivariate winsorization, a bivariate normal 
#' distribution is assumed and the data are shrunken towards the boundary of a 
#' tolerance ellipse with coverage probability \code{prob}.  The boundary of 
#' this ellipse is thereby given by all points that have a squared Mahalanobis 
#' distance equal to the quantile of the \eqn{\chi^{2}}{chi-squared} 
#' distribution given by \code{prob}.  Furthermore, the initial correlation 
#' matrix required for the Mahalanobis distances is computed based on adjusted 
#' univariate winsorization.
#' 
#' @param x  a numeric vector.
#' @param y  a numeric vector.
#' @param type  a character string specifying the type of winsorization to be 
#' used.  Possible values are \code{"univariate"} for univariate winsorization, 
#' \code{"adjusted"} for adjusted univariate winsorization, or 
#' \code{"bivariate"} for bivariate winsorization.
#' @param standardized  a logical indicating whether the data are already 
#' robustly standardized.
#' @param centerFun  a function to compute a robust estimate for the center to 
#' be used for robust standardization (defaults to 
#' \code{\link[stats]{median}}).  Ignored if \code{standardized} is \code{TRUE}.
#' @param scaleFun  a function to compute a robust estimate for the scale to 
#' be used for robust standardization (defaults to \code{\link[stats]{mad}}).  
#' Ignored if \code{standardized} is \code{TRUE}.
#' @param const  numeric; tuning constant to be used in univariate or adjusted 
#' univariate winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in bivariate
#' winsorization (defaults to 0.95).
#' @param tol  a small positive numeric value.  This is used in bivariate 
#' winsorization to determine whether the initial estimate from adjusted 
#' univariate winsorization is close to 1 in absolute value.  In this case, 
#' bivariate winsorization would fail since the points form almost a straight 
#' line, and the initial estimate is returned.
#' @param \dots  additional arguments to be passed to 
#' \code{\link[=standardize]{robStandardize}}.
#' 
#' @return The robust correlation estimate.
#' 
#' @author Andreas Alfons, based on code by Jafar A. Khan, Stefan Van Aelst and 
#' Ruben H. Zamar
#' 
#' @references 
#' Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust linear model 
#' selection based on least angle regression. \emph{Journal of the American 
#' Statistical Association}, \bold{102}(480), 1289--1299.
#' 
#' @seealso \code{\link{winsorize}}
#' 
#' @examples 
#' ## generate data
#' library("mvtnorm")
#' set.seed(1234)  # for reproducibility
#' Sigma <- matrix(c(1, 0.6, 0.6, 1), 2, 2)
#' xy <- rmvnorm(100, sigma=Sigma)
#' x <- xy[, 1]
#' y <- xy[, 2]
#' 
#' ## introduce outlier
#' x[1] <- x[1] * 10
#' y[1] <- y[1] * (-5)
#' 
#' ## compute correlation
#' cor(x, y)
#' corHuber(x, y)
#' 
#' @keywords multivariate robust
#' 
#' @export


# robust correlation based on winsorization
corHuber <- function(x, y, type = c("bivariate", "adjusted", "univariate"), 
                     standardized = FALSE, centerFun = median, 
                     scaleFun = mad, const = 2, prob = 0.95, 
                     tol = .Machine$double.eps^0.5, ...) {
  ## initializations
  n <- length(x)
  if(length(y) != n) stop("'x' and 'y' must have the same length")
  if(n == 0) return(NA)  # zero length vectors
  type <- match.arg(type)
  if(!isTRUE(standardized)) {
    # robustly standardize 'x'
    x <- robStandardize(x, centerFun=centerFun, scaleFun=scaleFun, ...)
    # robustly standardize 'y'
    y <- robStandardize(y, centerFun=centerFun, scaleFun=scaleFun, ...)
  }
  ## compute robust correlation according to winsorization type
  switch(type, 
         bivariate=corHuberBi(x, y, const=const, prob=prob, tol=tol),
         adjusted=corHuberAdj(x, y, const=const),
         univariate=corHuberUni(x, y, const=const))
}


## the following functions assume that the data are already standardized

# robust correlation based on bivariate winsorization
corHuberBi <- function(x, y, const = 2, prob = 0.95, 
                       tol = .Machine$double.eps^0.5) {
  # call C++ function
  .Call("R_corHuberBi", R_x=x, R_y=y, R_c=const, R_prob=prob, R_tol=tol, 
        PACKAGE="robustHD")
}

# robust correlation based on adjusted univariate winsorization
corHuberAdj <- function(x, y, const = 2) {
  # call C++ function
  .Call("R_corHuberAdj", R_x=x, R_y=y, R_c=const, PACKAGE="robustHD")
}

# robust correlation based on univariate winsorization
# with barebones Pearson correlation, C++ code is slightly faster than R code
corHuberUni <- function(x, y, const = 2) {
  # call C++ function
  .Call("R_corHuberUni", R_x=x, R_y=y, R_c=const, PACKAGE="robustHD")
}

# -------------------------------------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
#
# based on code by Jafar A. Khan, Stefan Van Aelst and Ruben H. Zamar
# -------------------------------------------------------------------

#' Data cleaning by winsorization
#' 
#' Clean data by means of winsorization, i.e., by shrinking outlying 
#' observations to the border of the main part of the data.
#' 
#' The borders of the main part of the data are defined on the scale of the 
#' robustly standardized data.  In the univariate case, the borders are given 
#' by \eqn{+/-}\code{const}, thus a symmetric distribution is assumed.  In the 
#' multivariate case, a normal distribution is assumed and the data are 
#' shrunken towards the boundary of a tolerance ellipse with coverage 
#' probability \code{prob}.  The boundary of this ellipse is thereby given by 
#' all points that have a squared Mahalanobis distance equal to the quantile of 
#' the \eqn{\chi^{2}}{chi-squared} distribution given by \code{prob}.
#' 
#' @param x  a numeric vector, matrix or data frame to be cleaned.
#' @param standardized  a logical indicating whether the data are already 
#' robustly standardized.
#' @param centerFun  a function to compute a robust estimate for the center to 
#' be used for robust standardization (defaults to 
#' \code{\link[stats]{median}}).  Ignored if \code{standardized} is \code{TRUE}.
#' @param scaleFun  a function to compute a robust estimate for the scale to 
#' be used for robust standardization (defaults to \code{\link[stats]{mad}}).  
#' Ignored if \code{standardized} is \code{TRUE}.
#' @param const  numeric; tuning constant to be used in univariate 
#' winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in multivariate
#' winsorization (defaults to 0.95).
#' @param tol  a small positive numeric value used to determine singularity 
#' issues in the computation of correlation estimates based on bivariate 
#' winsorization (see \code{\link{corHuber}}).
#' @param return  character string; if \code{standardized} is \code{TRUE}, 
#' this specifies the type of return value.  Possible values are \code{"data"} 
#' for returning the cleaned data, or \code{"weights"} for returning data 
#' cleaning weights.
#' @param \dots  for the generic function, additional arguments to be passed 
#' down to methods.  For the \code{"data.frame"} method, additional arguments 
#' to be passed down to the \code{"matrix"} method.  For the other methods, 
#' additional arguments to be passed down to 
#' \code{\link[=standardize]{robStandardize}}.
#' 
#' @return 
#' If \code{standardize} is \code{TRUE} and \code{return} is \code{"weights"}, 
#' a set of data cleaning weights.  Multiplying each observation of the 
#' standardized data by the corresponding weight yields the cleaned 
#' standardized data.
#' 
#' Otherwise an object of the same type as the original data \code{x} 
#' containing the cleaned data is returned.
#' 
#' @note Data cleaning weights are only meaningful for standardized data.  In 
#' the general case, the data need to be standardized first, then the data 
#' cleaning weights can be computed and applied to the standardized data, after 
#' which the cleaned standardized data need to be backtransformed to the 
#' original scale.
#' 
#' @author Andreas Alfons, based on code by Jafar A. Khan, Stefan Van Aelst and 
#' Ruben H. Zamar
#' 
#' @references 
#' Khan, J.A., Van Aelst, S. and Zamar, R.H. (2007) Robust linear model 
#' selection based on least angle regression. \emph{Journal of the American 
#' Statistical Association}, \bold{102}(480), 1289--1299.
#' 
#' @seealso \code{\link{corHuber}}
#' 
#' @examples 
#' ## generate data
#' set.seed(1234)     # for reproducibility
#' x <- rnorm(10)     # standard normal
#' x[1] <- x[1] * 10  # introduce outlier
#' 
#' ## winsorize data
#' x
#' winsorize(x)
#' 
#' @keywords robust
#' 
#' @export

winsorize <- function(x, ...) UseMethod("winsorize")


#' @rdname winsorize
#' @method winsorize default
#' @export

winsorize.default <- function(x, standardized = FALSE, centerFun = median, 
                              scaleFun = mad, const = 2, 
                              return = c("data", "weights"), ...) {
  ## initializations
  standardized <- isTRUE(standardized)
  if(standardized) return <- match.arg(return)
  else {
    # standardize data
    x <- robStandardize(x, centerFun=centerFun, scaleFun=scaleFun, ...)
    center <- attr(x, "center")
    scale <- attr(x, "scale")
  }
  ## winsorize standardized data
#   ind <- abs(x) > const           # observations in 'x' that need to be shrunken
#   x[ind] <- const * sign(x[ind])  # winsorize
  weights <- pmin(const / abs(x), 1)
  if(standardized && return == "weights") return(weights)
  x <- weights * x
  ## finalizations
  if(!standardized) {
    # transform back to original scale and remove attributes
    x <- c(x * scale + center)
  }
  x
}


#' @rdname winsorize
#' @method winsorize matrix
#' @export

winsorize.matrix <- function(x, standardized = FALSE, centerFun = median, 
                             scaleFun = mad, const = 2, prob = 0.95, 
                             tol = .Machine$double.eps^0.5, 
                             return = c("data", "weights"), ...) {
  ## initializations
  m <- ncol(x)
  standardized <- isTRUE(standardized)
  if(standardized) return <- match.arg(return)
  ## winsorize the data
  if(m == 1) {
    attributes <- attributes(x)
    x <- winsorize(c(x), standardized=standardized, centerFun=centerFun, 
                   scaleFun=scaleFun, const=const, return=return, ...)
    if(!(standardized && return == "weights")) attributes(x) <- attributes
  } else {
    if(nrow(x) <= m) {
      stop("not enough observations for inversion of correlation matrix")
    }
    if(!standardized) {
      # standardize data
      attributes <- attributes(x)
      x <- robStandardize(x, centerFun=centerFun, scaleFun=scaleFun, ...)
      center <- attr(x, "center")
      scale <- attr(x, "scale")
    }
    ## compute correlation matrix (call C++ function)
    R <- .Call("R_corMatHuber", R_x=x, R_c=const, R_prob=prob, R_tol=tol, 
               PACKAGE="robustHD")
    # check if correlation matrix is positive definite and perform 
    # corrections if necessary
    eig <- eigen(R, symmetric=TRUE)
    if(eig$values[m] < 0) {  # last eigenvalue is the smallest
      Q <- eig$vectors  # matrix of eigenvectors
      lambda <- apply(x %*% Q, 2, scaleFun)^2
      R <- Q %*% diag(lambda) %*% t(Q)
    }
    ## compute Mahalanobis distances
    # function 'mahalanobis()' requires suppling estimates for the centers, 
    # hence it is not used here
    invR <- solve(R)
    md <- rowSums((x %*% invR) * x)  # squared Mahalanobis distance
    ## shrink observations with too large distances
    d <- qchisq(prob, df=m)  # quantile of the chi-squared distribution
    #        ind <- which(md > d)  # select observations that need to be shrunken
    #        x[ind,] <- x[ind,] * sqrt(d / md[ind])  # shrink selected observations
    weights <- pmin(sqrt(d / md), 1)
    if(standardized && return == "weights") return(weights)
    x <- weights * x
    ## finalizations
    if(!standardized) {
      # transform back to original scale and remove attributes
      x <- sweep(x, 2, scale, "*", check.margin=FALSE)
      x <- sweep(x, 2, center, "+", check.margin=FALSE)
      attributes(x) <- attributes
    }
  }
  x
}


#' @rdname winsorize
#' @method winsorize data.frame
#' @export

winsorize.data.frame <- function(x, ...) winsorize(as.matrix(x, ...))

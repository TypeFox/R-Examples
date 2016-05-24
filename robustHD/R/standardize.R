# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Data standardization
#' 
#' Standardize data with given functions for computing center and scale.
#' 
#' @param x  a numeric vector, matrix or data frame to be standardized.
#' @param centerFun  a function to compute an estimate of the center of a 
#' variable (defaults to \code{\link{mean}}).
#' @param scaleFun  a function to compute an estimate of the scale of a 
#' variable (defaults to \code{\link[stats]{sd}}).
#' 
#' @return An object of the same type as the original data \code{x} containing 
#' the centered and scaled data.  The center and scale estimates of the 
#' original data are returned as attributes \code{"center"} and \code{"scale"}, 
#' respectively.
#' 
#' @note The implementation contains special cases for the typically used 
#' combinations \code{\link{mean}}/\code{\link[stats]{sd}} and 
#' \code{\link[stats]{median}}/\code{\link[stats]{mad}} in order to reduce 
#' computation time.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{scale}}, \code{\link{sweep}}
#' 
#' @examples 
#' ## generate data
#' set.seed(1234)     # for reproducibility
#' x <- rnorm(10)     # standard normal
#' x[1] <- x[1] * 10  # introduce outlier
#' 
#' ## standardize data
#' x
#' standardize(x)     # mean and sd
#' robStandardize(x)  # median and MAD
#' 
#' @keywords array
#' 
#' @export

standardize <- function(x, centerFun = mean, scaleFun = sd) {
  if(is.null(dim(x))) {
    # generic standardization with special cases mean/sd and median/MAD
    center <- centerFun(x)  # compute center
    x <- x - center  # sweep out center
    if(identical(centerFun, mean) && identical(scaleFun, sd)) {
      # classical standardization with mean and standard deviation
      scale <- sqrt(sum(x^2) / max(1, length(x)-1))
    } else if(identical(centerFun, median) && identical(scaleFun, mad)) {
      # robust standardization with median and MAD
      # compute MAD with median already swept out
      scale <- mad(x, center=0)
    } else scale <- scaleFun(x)  # compute scale
    x <- x / scale  # sweep out scale
  } else {
    # generic standardization with special cases mean/sd and median/MAD
    if(identical(centerFun, mean)) {
      center <- colMeans(x)  # compute column means (faster than apply)
    } else center <- apply(x, 2, centerFun)  # compute column centers
    x <- sweep(x, 2, center, check.margin=FALSE)  # sweep out column centers
    if(identical(centerFun, mean) && identical(scaleFun, sd)) {
      # classical standardization with mean and standard deviation
      f <- function(v) sqrt(sum(v^2) / max(1, length(v)-1))
      scale <- apply(x, 2, f)
    } else if(identical(centerFun, median) && identical(scaleFun, mad)) {
      # robust standardization with median and MAD
      # compute column MADs with median already swept out
      scale <- apply(x, 2, mad, center=0)
    } else scale <- apply(x, 2, scaleFun)  # compute column scales
    x <- sweep(x, 2, scale, "/", check.margin=FALSE)  # sweep out column scales
  }
  # add attributes and return standardized data
  attr(x, "center") <- center
  attr(x, "scale") <- scale
  x
}


# wrapper function to robustly standardize data
# if the robust scale of a variable is too small, it is standardized with mean 
# and standard deviation

#' @rdname standardize
#' 
#' @param fallback  a logical indicating whether standardization with 
#' \code{\link{mean}} and \code{\link[stats]{sd}} should be performed as a 
#' fallback mode for variables whose robust scale estimate is too small.  This 
#' is useful, e.g., for data containing dummy variables.
#' @param eps  a small positive numeric value used to determine whether the 
#' robust scale estimate of a variable is too small (an effective zero).
#' @param \dots  currently ignored.
#' 
#' @details 
#' \code{robStandardize} is a wrapper function for robust standardization, 
#' hence the default is to use \code{\link[stats]{median}} and 
#' \code{\link[stats]{mad}}.
#' 
#' @export

robStandardize <- function(x, centerFun = median, scaleFun = mad, 
                           fallback = FALSE, eps = .Machine$double.eps, ...) {
  # robustly standardize data
  xs <- standardize(x, centerFun=centerFun, scaleFun=scaleFun)
  # if requested, check if some variables have too small a robust scale and 
  # use standardization with mean/sd as fallback mode
  if(isTRUE(fallback)) {
    scale <- attr(xs, "scale")
    if(is.null(dim(x))) {
      if(scale <= eps) xs <- standardize(x)
    } else {
      tooSmall <- which(scale <= eps)
      if(length(tooSmall) > 0) {
        # standardize with mean and standard deviation
        center <- attr(xs, "center")
        xcs <- standardize(x[, tooSmall])
        center[tooSmall] <- attr(xcs, "center")
        scale[tooSmall] <- attr(xcs, "scale")
        xs[, tooSmall] <- xcs
        attr(xs, "center") <- center
        attr(xs, "scale") <- scale
      }
    }
  }
  # return standardized data
  xs
}

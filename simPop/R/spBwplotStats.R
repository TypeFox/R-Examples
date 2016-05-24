#' Weighted box plot statistics
#' 
#' Compute the statistics necessary for producing box-and-whisker plots of
#' continuous or semi-continuous variables, taking into account sample weights.
#' 
#' The function \code{\link{quantileWt}} is used for the computation of
#' (weighted) quantiles.  The median is computed together with the first and
#' the third quartile, which form the box.  If \code{range} is positive, the
#' whiskers extend to the most extreme data points that have a distance to the
#' box of no more than \code{coef} times the interquartile range.  For
#' \code{coef = 0}, the whiskers mark the minimum and the maximum of the
#' sample, whereas a negative value causes an error.
#' 
#' @name spBwplotStats
#' @param x a numeric vector.
#' @param weights an optional numeric vector containing sample weights.
#' @param coef a numeric value that determines the extension of the whiskers.
#' @param zeros a logical indicating whether the variable specified by
#' \code{additional} is semi-continuous, i.e., contains a considerable amount
#' of zeros.  If \code{TRUE}, the (weighted) box plot statistics are computed
#' for the non-zero data points only and the number of zeros is returned, too.
#' @param do.out a logical indicating whether data points that lie beyond the
#' extremes of the whiskers should be returned.
#' @return A list of class \code{"spBwplotStats"} with the following
#' components: \item{stats}{A vector of length 5 containing the (weighted)
#' statistics for the construction of a box plot.} \item{n}{if \code{weights}
#' is \code{NULL}, the number of non-missing and, if \code{zeros} is
#' \code{TRUE}, non-zero data points.  Otherwise the sum of the weights of the
#' corresponding points.} \item{nzero}{if \code{zeros} is \code{TRUE} and
#' \code{weights} is \code{NULL}, the number of zeros.  If \code{zeros} is
#' \code{TRUE} and \code{weights} is not \code{NULL}, the sum of the weights of
#' the zeros.  If \code{zeros} is not \code{TRUE}, this is \code{NULL}.}
#' \item{out}{if \code{do.out}, the values of any data points that lie beyond
#' the extremes of the whiskers.}
#' @author Stefan Kraft and Andreas Alfons
#' @export
#' @seealso \code{\link{spBwplot}}, for producing (weighted) box plots of
#' continuous or semi-continuous variables.
#' 
#' \code{\link{quantileWt}} for the computation of (weighted) sample quantiles.
#' 
#' \code{\link[grDevices]{boxplot.stats}} for the unweighted statistics for box
#' plots (not considering semi-continuous variables).
#' @keywords dplot
#' @examples
#' 
#' data(eusilcS)
#' 
#' ## semi-continuous variable
#' spBwplotStats(eusilcS$netIncome, 
#'     weights=eusilcS$rb050, do.out = FALSE)
#' 
spBwplotStats <- function(x, weights = NULL, coef = 1.5, zeros = TRUE, do.out = TRUE) {
  # initializations
  if ( !is.numeric(x) ) {
    stop("'x' must be a numeric vector!\n")
  }
  if ( !is.numeric(coef) || length(coef) != 1 || coef < 0 ) {
    stop("'coef' must be a single non-negative number!\n")
  }
  # get quantiles
  if ( isTRUE(zeros) ) {
    zero <- ifelse(is.na(x), FALSE, x == 0)
    x <- x[!zero]
    if ( is.null(weights) ) {
      nzero <- sum(zero)
    } else {
      # if 'zeros' is not TRUE, these checks are done in 'quantileWt'
      # but here we need them since we use subscripting
      if ( !is.numeric(weights) ) {
        stop("'weights' must be a numeric vector!\n")
      }
      else if ( length(weights) != length(zero) ) {
        stop("'weights' must have the same length as 'x'!\n")
      }
      nzero <- sum(weights[zero])
      weights <- weights[!zero]
    }
  } else {
    nzero <- NULL
  }
  ok <- !is.na(x)
  n <- if ( is.null(weights) ) {
    sum(ok)
  } else {
    sum(weights[ok])
  }
  if ( n == 0 ) {
    stats <- rep.int(NA, 5)
  } else {
    stats <- quantileWt(x, weights)
  }
  iqr <- diff(stats[c(2, 4)])  # inter quartile range
  if ( coef == 0 ) {
    do.out <- FALSE
  } else {
    if ( is.na(iqr) ) {
      out <- is.infinite(x)
    } else {
      lower <- stats[2] - coef * iqr
      upper <- stats[4] + coef * iqr
      out <- ifelse(ok, x < lower | x > upper, FALSE)
    }
    if ( any(out) ) {
      stats[c(1, 5)] <- range(x[!out], na.rm=TRUE)
    }
  }
  res <- list(stats=stats, n=n, nzero=nzero, out=if(isTRUE(do.out)) x[out] else numeric())
  class(res) <- "spBwplotStats"
  res
}


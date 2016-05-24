# Created by Matthew A. Birk
# Calculates the standard error of a vector
# Last updated: Jun 2014

#' Standard Error
#'
#' Computes the standard error of the values in \code{x}. If \code{na.rm} is TRUE then missing values are removed before computation proceeds.
#'
#' @param x a numeric vector or an R object which is coercible to one by as.vector(x, "numeric”).
#' @param na.rm logical. Should missing values be removed?
#'
#' @author Matthew A. Birk, \email{matthewabirk@@gmail.com}
#' @seealso \code{\link{sd}}, \code{\link{var}}
#'
#' @examples
#' se(1:10)
#'
#' @encoding UTF-8
#' @export
#' @import stats

se=function(x,na.rm=FALSE)
{
  return(sqrt(stats::var(if (is.vector(x)) x else as.double(x), na.rm = na.rm) / length(x)))
}
#' Sample range
#'
#' This function computes the sample range.
#'
#' "The range is the difference between the largest number and the smallest
#' number in the set." Source: Onwubiko page 176.
#'
#' The following statements are from \code{\link[base]{range}}:
#'
#' "If na.rm is \code{FALSE}, \code{NA} and \code{NaN} values in any of the
#' arguments will cause \code{NA} values to be returned, otherwise \code{NA}
#' values are ignored."
#'
#' "If finite is \code{TRUE}, the minimum and maximum of all finite values is
#' computed, i.e., \code{finite = TRUE} includes \code{na.rm = TRUE}."
#'
#'
#'
#' @param x numeric vector
#' @param na.rm logical vector that determines whether the missing values
#'    should be removed or not.
#' @param finite logical vector that determines whether non-finite values
#'    should be removed or not.
#'
#' @return ranges as the difference between the maximum and minimum values in \code{x}
#'    as a numeric \code{\link{vector}}. Unlike the \code{\link[base]{range}}, ranges can't
#'    take character vectors as arguments, only numeric vectors.
#'
#'
#' @references
#' Chinyere Onwubiko, \emph{An Introduction to Engineering}, Mission, Kansas: Schroff Development Corporation, 1997, page 176.
#'
#'
#' @encoding UTF-8
#'
#' @seealso \code{\link{sgm}} for geometric mean, \code{\link{shm}} for harmonic mean, \code{\link{cv}} for
#'  coefficient of variation (CV), \code{\link{rms}} for root-mean-square (RMS), \code{\link{relerror}}
#'  for relative error, and \code{\link{approxerror}} for approximate error.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @examples
#' library(iemisc)
#' require(stats)
#' set.seed(100) # makes the example reproducible
#' x <- rnorm(100)
#' ranges(x)
#'
#'
#'
#' @export
ranges <- function (x, na.rm = FALSE, finite = FALSE) {

# The moments::kurtosis code has been helpful with regards to the treatment of na.rm

newx <- range(x, na.rm = na.rm, finite = finite)
# R's range (minimum and maximum value)

ranges <- diff(newx)
# sample range

return(ranges)
}

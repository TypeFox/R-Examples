#' Constrain Data to a Specified Range
#'
#' Constrain data to a specified range, assigning values from the specified
#'   range to those outside the range, typically for graphing purposes.
#' @param x
#'   A numeric vector of values to constrain.
#' @param xrange
#'   A numeric vector of length two specifying the constraints, the minimum and
#'     maximum value for \code{x}.
#' @return
#'   A numeric vector, the same length as \code{x}, in which the minimum
#'     constraint is assigned to values of \code{x} less than the minimum,
#'     and the maximum constraint is assigned to values of \code{x} greater than
#'     the maximum.
#' @export
#' @examples
#' constrain(1:20, c(3, 19))

constrain <- function(x, xrange) {
  if (!is.numeric(x)) stop("x must be numeric.")
  if (length(xrange)!=2 | any(is.na(xrange)) | !is.numeric(xrange)) {
    stop("xrange must be a non-missing numeric vector of length 2")
  }
  x[!is.na(x) & x<xrange[1]] <- xrange[1]
  x[!is.na(x) & x>xrange[2]] <- xrange[2]
  x
  }

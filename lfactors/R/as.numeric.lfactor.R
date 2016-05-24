#' Numeric Vectors from lfactors
#' @method as.numeric lfactor
#' @description
#' Returns numeric representation of an lfactor equal to the \code{levels} argument for each value. This is
#' different from the behavior of factor which would ignore the values of \code{level}.
#' 
#' @param x same as \code{\link[base]{as.numeric}}
#' @param \dots not used
#' @details
#' This method does not return floating point (numeric) results that are otherwise equal to the results from \code{\link{as.integer.lfactor}}.
#' Instead it returns the value of the level that was input when the lfactor was created.
#'
#' @seealso \code{\link[base]{as.numeric}}, \code{\link{as.integer.lfactor}}
#' @example man/examples/as.num.R
#' @export
as.numeric.lfactor <- function(x, ...) {
  as.numeric(as.character(switchllevels(x)))
}


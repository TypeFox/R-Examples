#' Integer Vectors from lfactors
#' @method as.integer lfactor
#' @description
#' Returns integer representation of an lfactor that ignores the values used in the \code{levels} argument
#' when the lfactor was created and instead returns an integer representation starting with 1.
#' 
#' @param x same as \code{\link[base]{as.integer}}
#' @param \dots not used
#' @details
#' This method does not return integer results that are otherwise equal to the results from as.numeric
#' for compatability with \code{\link[Matrix]{sparse.model.matrix}}.
#'
#' @seealso \code{\link[base]{as.integer}}, \code{\link{as.numeric.lfactor}}
#' @example man/examples/as.num.R
#' @export
as.integer.lfactor <- function(x, ...) {
  as.integer(as.factor.lfactor(x))
  #as.integer(as.character(switchllevels(x)))
}

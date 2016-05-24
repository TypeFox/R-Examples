#' factor from an lfactor
#' @method as.factor lfactor
#' @description
#' Returns an \code{\link[base]{factor}} from an \code{\link{lfactor}}.
#' 
#' @param x the lfactor to be coerced to a factor
#' @details
#' simply drops the numeric levels from the lfactor and returns a normal factor.
#'
#' @seealso \code{\link[base]{as.factor}}
#' @export
as.factor.lfactor <- function(x) {
  class(x) <- "factor"
  attr(x, "llevels") <- NULL
  x
}

methods::setGeneric("as.factor")
methods::setMethod("as.factor", methods::signature(x="lfactor"), as.factor.lfactor)

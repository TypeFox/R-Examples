#' @title create method xtree.pls
#' @details
#' Internal function. 
#' @param x the element representing the method.
#' @param \dots Further arguments passed on to \code{\link{xtree.pls}}.
#' @return internal method
#' @keywords internal
#' @export

xtree.pls <- function(x, ...) UseMethod("xtree.pls")

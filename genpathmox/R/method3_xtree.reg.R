#' @title create method xtree.reg
#' @details
#' Internal function. 
#' @param x the element representing the method.
#' @param \dots Further arguments passed on to \code{\link{xtree.reg}}.
#' @return internal method
#' @keywords internal
#' @export

xtree.reg <- function(x, ...) UseMethod("xtree.reg")

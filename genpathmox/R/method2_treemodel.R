#' @title create method treemodel 
#' @details
#' Internal function. 
#' @param x the element representing the method.
#' @param \dots Further arguments passed on to \code{\link{treemodel}}.
#' @return internal method
#' @keywords internal
#' @export

treemodel <- function(x, ...) UseMethod("treemodel")

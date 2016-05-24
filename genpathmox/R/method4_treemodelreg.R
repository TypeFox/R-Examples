#' @title create method treemodelreg
#' @details
#' Internal function. 
#' @param x the element representing the method.
#' @param \dots Further arguments passed on method.
#' @return internal method
#' @keywords internal
#' @export

treemodelreg <- function(x, ...) UseMethod("treemodelreg")

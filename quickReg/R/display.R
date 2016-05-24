#' Generic function to display a data.frame or `reg` object
#'

#' \code{display} is a generic function to display summary information of varibles in a data.frame or summary result of `reg` class
#' @param x an object to display, a data.frame or a  `reg` class
#' @param \dots additional arguments

#' @export
#' @seealso \code{\link{display.data.frame}}, \code{\link{display.reg}}

display <- function(x, ...) UseMethod("display") 

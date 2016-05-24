##' Generic function \code{color}
##'
##' This generic function supports extraction of colors from colored
##' objects
##' @param x A colored object
##' @param \dots Additional arguments
##' @return A vector of colors
##' @rdname color
##' @export color
color <- function(x, ...) UseMethod("color")

##' Generic function 'key'
##'
##' This generic function supports extraction of color keys from colored
##' factors, via the specific method \code{key.colored}
##' @param x A colored object
##' @param \dots Additional arguments
##' @return A vector of colors
##' @rdname key
##' @export
key <- function(x, ...) UseMethod("key")

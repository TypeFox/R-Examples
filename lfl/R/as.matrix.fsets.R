as.matrix.fsets <- function(x, ...) {
    if (!is.fsets(x)) {
        stop("'x' must be an instance of the 'fsets' class")
    }
    class(x) <- 'matrix'
    attr(x, 'vars') <- NULL
    attr(x, 'specs') <- NULL
    return(x)
}

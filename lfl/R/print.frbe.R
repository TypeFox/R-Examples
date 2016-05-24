print.frbe <- function(x, ...) {
    if (!is.frbe(x)) {
        stop("'x' must be an instance of the 'frbe' class")
    }
    xx <- x
    class(xx) <- NULL
    print(xx, ...)
}


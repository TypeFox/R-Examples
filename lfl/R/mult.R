mult <- function(x, y, f, ...) {
    stopifnot(is.matrix(x))
    stopifnot(is.matrix(y))
    stopifnot(nrow(x) > 0)
    stopifnot(ncol(y) > 0)
    stopifnot(ncol(x) == nrow(y))
    stopifnot(is.function(f))
    stopifnot(length(formals(f)) >= 2) # f should have at least 2 arguments 

    ff <- function(xx, yy) {
        f(xx, yy, ...)
    }
    .Call('mult', x, y, ff, PACKAGE='lfl')
}

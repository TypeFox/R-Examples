qmatmult <- function(x, y) {

    if (! (is.numeric(x) || is.character(x)))
        stop("'x' not numeric or character")
    if (! (is.numeric(y) || is.character(y)))
        stop("'y' not numeric or character")

    if (is.numeric(x))
        x <- d2q(x)
    if (is.numeric(y))
        y <- d2q(y)

    if (! (is.character(x) && is.character(y))) {
        stop("Cannot happen!")
    }

    if (! (is.matrix(x) && is.matrix(y)))
        stop("'x' and 'y' must be matrices")

    stopifnot(ncol(x) == nrow(y))

    .Call("qmatmult", x, y, PACKAGE = "rcdd")
}

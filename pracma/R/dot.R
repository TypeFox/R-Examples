###
### DOT.R  Scalar product
###

dot <- function(x, y) {
    if (length(x) == 0 && length(y) == 0) return(0)
    if (!(is.numeric(x) || is.complex(x)) ||
        !(is.numeric(y) || is.complex(y)))
        stop("Arguments 'x' and 'y' must be real or complex.")
    x <- drop(x); y <- drop(y)
    if (any(dim(x) != dim(y)))
        stop("Matrices 'x' and 'y' must be of same size")

    if (is.vector(x) && is.vector(y)) {
        dim(x) <- c(length(x), 1)
        dim(y) <- c(length(y), 1)
    }
    x.y <- apply(Conj(x) * y, 2, sum)
    return(x.y)
}

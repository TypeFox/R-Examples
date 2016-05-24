##
##  e p s . R
##


eps <- function(x = 1.0) {
    stopifnot(is.numeric(x))

    x <- max(abs(x))

    if (x <  .Machine$double.xmin) {
        e <- .Machine$double.xmin
    } else {
        e <- 2^floor(log2(x)) * .Machine$double.eps
    }
    e
}

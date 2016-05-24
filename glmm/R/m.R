
matvecmult <- function(a, x)
{
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(x))
    stopifnot(is.finite(a))
    stopifnot(is.finite(x))
    stopifnot(is.matrix(a))
    stopifnot(ncol(a) == length(x))
    .C("matvecmult", a = as.double(a), x = as.double(x),
        nrow = nrow(a), ncol = ncol(a), result = double(nrow(a)),
        PACKAGE = "glmm")$result
}

matmatmult <- function(a, b)
{
    stopifnot(is.numeric(a))
    stopifnot(is.numeric(b))
    stopifnot(is.finite(a))
    stopifnot(is.finite(b))
    stopifnot(is.matrix(a))
    stopifnot(is.matrix(b))
    stopifnot(ncol(a) == nrow(b))
    .C("matmatmult", a = as.double(a), b = as.double(b),
        nrowa = nrow(a), ncola = ncol(a), ncolb = ncol(b),
        c = matrix(as.double(0), nrow = nrow(a), ncol = ncol(b)),
        PACKAGE = "glmm")$c
}


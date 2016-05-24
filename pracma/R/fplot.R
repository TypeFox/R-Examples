##
##  f p l o t . R
##


fplot <- function(f, a, b, n = 101, vectorized = FALSE,
                  col = 1:8, lty = 1, lwd = 1, ...) {
    # f is a multivariate function, vectorized in columns
    x <- linspace(a, b, n)
    m <- length(f(a))
    if (length(f(b)) != m)
        stop("Function 'f' is not truly multivariate: length(f(a)) != length(f(b)).")
    if (size(f(c(a,b)), 1) != 2 && vectorized) {
        warning("Function 'f' not correctly vectorized: size(f(c(a,b)), 1) != 2.")
        vectorized <- FALSE
    }

    if (!vectorized) {
        A <- matrix(NA, n, m)
        for (i in 1:n) {
            A[i, ] <- f(x[i])
        }
    } else {
        A <- f(x)
    }

    matplot(A, type = "l", col = col, lty = lty, lwd = lwd, ...)
    grid()
    invisible()
}

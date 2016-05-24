##
##  a i t k e n . R  Aitken's acceleration method
##


aitken <- function(f, x0, nmax = 12, tol = 1e-8, ...) {
    if (!is.numeric(x0) || length(x0) != 1)
        stop("Argument 'x0' must be a numeric scalar.")
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    x <- x0
    diff <- 1 + tol
    niter <- 0
    while (diff > tol && niter <= nmax) {
        gx <- f(x)
        ggx <- f(gx)
        xnew <- (x*ggx - gx^2) / (ggx - 2*gx + x)
        diff <- abs(x - xnew)
        x <- xnew
        niter <- niter + 1
    }
    if (niter > nmax)
        warning("Maximum number of iterations exceeded.")
    return(x)
}

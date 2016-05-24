##
##  l a g u e r r e . R  Laguerre Method
##


laguerre <- function(p, x0, nmax = 25, tol = .Machine$double.eps^(1/2)) {
    if (!is.numeric(p) && !is.complex(p))
        stop("Argument 'p' must be a numeric or complex vector.")
    if ( (!is.numeric(x0) && !is.complex(x0)) || length(x0) != 1)
        stop("Argument 'x0' must be a real or complex number.")

    n <- length(p) - 1
    p1 <- polyder(p)
    p2 <- polyder(p1)

    y0 <- polyval(p, x0)
    if (abs(y0) < tol) return(x0)

    for (m in 1:nmax) {
        a <- polyval(p1, x0) / y0
        a2 <- a^2
        b <- a2 - polyval(p2, x0) / y0

        x <- x0 - n/(a + a/abs(a) * sqrt((n-1)*(n*b - a2)))
        y <- polyval(p, x)
        if (abs(y) < tol) break
        x0 <- x
        y0 <- y
    }
    if (m > nmax)
        warning("Root finding process might not have converged.")

    return(x)
}

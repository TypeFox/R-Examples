##
##  s t e e p _ d e s c e n t . R  Minimization by Steepest Descent
##


steep_descent <- function (x0, f, g = NULL, info = FALSE,
                           maxiter = 100, tol = .Machine$double.eps^(1/2)) {
    eps <- .Machine$double.eps
    if (! is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.")
    n <- length(x0)

    # User provided or numerical gradient
    f <- match.fun(f)
    if (is.null(g)) g <- function(x) grad(f, x)
    else            g <- match.fun(g)

    if (info) cat(0, "\t", x0, "\n")

    x <- x0
    k <- 1
    while (k <= maxiter) {
        f1 <- f(x)
        g1 <- g(x)
        z1 <- sqrt(sum(g1^2))
        if (z1 == 0) {
            warning(
                paste("Zero gradient at:", x, f1, "-- not applicable.\n"))
            return(list(xmin = NA, fmin = NA, niter = k))
        }
        # else use gradient as unit vector
        g1 <- g1 / z1

        a1 <- 0
        a3 <- 1; f3 <- f(x - a3*g1)

        # Find a minimum on the gradient line (or line search)
        while (f3 >= f1) {
            a3 <- a3/2; f3 <- f(x - a3*g1)
            if (a3 < tol/2) {
                if (info)
                    cat("Method of steepest descent converged to:", x, "\n")
                x[x < eps] <- 0
                return(list(xmin = x, fmin = f(x), niter = k))
            }
        }

        # Check an intermediate point (for faster convergence)
        a2 <- a3/2; f2 <- f(x - a2*g1)
        h1 <- (f2 - f1)/a2
        h2 <- (f3 -f2)/(a3 - a2)
        h3 <- (h2 - h1)/a3
        a0 <- 0.5*(a2 - h1/h3); f0 <- f(x - a0*g1)

        if (f0 < f3) a <- a0
        else         a <- a3

        x <- x - a*g1
        if (info) cat(k, "\t", x, "\n")
        k <- k + 1
    }
    if(k > maxiter)
        warning("Maximum number of iterations reached -- not converged.\n")
    return(list(xmin = NA, fmin = NA, niter = k))
}

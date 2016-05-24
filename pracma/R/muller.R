##
##  m u l l e r . R  Muller's Method
##


muller <- function(f, p0, p1, p2 = NULL, maxiter = 100, tol = 1e-10) {
    if (is.null(p2)) p2 <- (p0 + p1)/2
    stopifnot(is.numeric(p0) || is.complex(p0), length(p0) == 1,
              is.numeric(p1) || is.complex(p1), length(p1) == 1,
              is.numeric(p2) || is.complex(p2), length(p2) == 1)
    f <- match.fun(f)

    fp0 <- f(p0); fp1 <- f(p1); fp2 <- f(p2)
    if (!is.finite(fp0) || !is.finite(fp1) || !is.finite(fp2))
        stop("Function 'f' not finite at one of the initial points.")

    # Initialization
    h1 <- p1 - p0
    h2 <- p2 - p1
    d1 <- (fp1 - fp0) / h1
    d2 <- (fp2 - fp1) / h2
    d  <- (d2 - d1) / (h2 + h1)

    # main loop
    i  <- 3
    while (i <= maxiter) {
        b <- d2 + h2*d
        D <- sqrt(b^2 - 4*f(p2)*d + 0i)

        if (abs(b - D) < abs(b + D)) {
        	E <- b + D
        } else {
            E <- b - D
        }

        h <- -2*f(p2) / E
        p <- p2 + h
        fp <- f(p)

        if (abs(h) < tol) break

        # prepare for next iteration
        p0 <- p1
        p1 <- p2
        p2 <- p
        h1 <- p1 - p0
        h2 <- p2 - p1
        d1 <- (f(p1) - f(p0)) / h1
        d2 <- (f(p2) - f(p1)) / h2
        d  <- (d2 - d1) / (h2 + h1)
        i  <- i + 1
    }

    if (i > maxiter)
        warning("Root not found to the desired tolerance.")

    if (abs(Im(p)) <= 0.1*tol ) {
        p <- Re(p); fp <- f(p)
    }

    return(list(root = p, fval = fp, niter = i, reltol = abs(h)))
}

##
##  s i m p a d p t . R  Adaptive Simpson's Rule
##


simpadpt <- function(f, a, b, tol = 1e-6, ...) {
    stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1)
    eps <- .Machine$double.eps

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)
    
    if (a == b)     return(0)
    else if (a > b) return(-1 * simpadpt(f, b, a, tol = tol))
        
    # Start with unequal subintervals
    h <- 1/8 * (b-a)
    x <- c(a,      a+h,    a+2*h,    (a+b)/2,    b-2*h,    b-h,    b)
    y <- c(f(a), f(a+h), f(a+2*h), f((a+b)/2), f(b-2*h), f(b-h), f(b))

    # Avoid infinities at end points
    if ( !is.finite(y[1]) ) y[1] <- f(a + eps*(b-a))
    if ( !is.finite(y[7]) ) y[7] <- f(b - eps*(b-a))

    # Call the adaptive simpson function
    hmin <- eps * (b-a) / 1024
    Q1 <- .simpadpt(f, x[1], x[3], y[1], y[2], y[3], tol, hmin)
    Q2 <- .simpadpt(f, x[3], x[5], y[3], y[4], y[5], tol, hmin)
    Q3 <- .simpadpt(f, x[5], x[7], y[5], y[6], y[7], tol, hmin)

    return(Q1 + Q2 + Q3)
}


.simpadpt <- function(f, a, b, fa, fc, fb, tol, hmin) {
    h <- b - a
    g <- (a + b)/2  # fc
    d <- (a + g)/2; fd <- f(d)
    e <- (g + b)/2; fe <- f(e)

    # Three- and five-point Simpson's rule
    # plus a one-step Romberg extrapolation
    Q1 <- (h/6) * (fa + 4*fc + fb)
    Q2 <- (h/12) * (fa + 4*fd + 2*fc + 4*fe + fb)
    Q <- Q2 + (Q2 - Q1)/15

    if (!is.finite(Q)) {
        warning("Infinite or NA function value encountered.")
        return(Q)
    } else if (abs(Q2 - Q) <= tol) {
        return(Q)
    } else if (abs(h) < hmin || g == a || g == b) {
        warning("Minimum step size reached; singularity possible.")
        return(Q)
    }

    Q4 <- .simpadpt(f, a, g, fa, fd, fc, tol, hmin)
    Q5 <- .simpadpt(f, g, b, fc, fe, fb, tol, hmin)

    return(Q4 + Q5)
}

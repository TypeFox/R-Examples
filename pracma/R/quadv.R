quadv <- function(f, a, b, tol = .Machine$double.eps^(1/2), ...) {
    stopifnot(is.numeric(a), is.numeric(b))
    if (length(a) != 1 || length(b) != 1 || a > b)
        stop("Interval boundaries must satisfy 'a <= b'.")
    if (a == b) return(list(Q = 0, fcnt = 0, estim.prec = 0))
    eps <- .Machine$double.eps

    fun <- match.fun(f)
    f   <- function(x) fun(x, ...)

    # Initialize with unequal intervals
    h <- 0.13579 * (b-a)
    x <- c(a, a + h, a + 2*h, (a + b)/2, b - 2*h, b - h, b)
    y <- f(x[1])
    for (j in 2:7) {
        y <- rbind(y, f(x[j]))
    }
    fcnt <- 7

    # Fudge endpoints to avoid infinities
    if (any(!is.finite(y[1, ]))) {
        y[1, ] <- f(a + eps*(b-a))
        fcnt <- fcnt + 1
    }
    if (any(!is.finite(y[7, ]))) {
        y[7, ] <- f(b - eps*(b-a))
        fcnt <- fcnt + 1
    }

    # Call recursively the main integrator function
    hmin <- eps*(b-a)/1024
    I1 <- .quadvstep(f,x[1], x[3], y[1, ],y[2, ], y[3, ], tol, fcnt, hmin)
    Q1 <- I1$Q
    fcnt <- I1$fcnt
    I2 <- .quadvstep(f,x[3], x[5], y[3, ],y[4, ], y[5, ], tol, fcnt, hmin)
    Q2 <- I2$Q
    fcnt <- I2$fcnt
    I3 <- .quadvstep(f,x[5], x[7], y[5, ],y[6, ], y[7, ], tol, fcnt, hmin)
    Q3 <- I3$Q
    fcnt <- I3$fcnt

    Q <- unname(Q1 + Q2 + Q3)
    return(list(Q = Q, fcnt = fcnt, estim.prec = tol*(fcnt-7)/2))
}


.quadvstep <- function(f, a, b, fa, fc, fb, tol, fcnt, hmin) {
    maxfcnt <- 10000

    # Evaluate integrand twice in interior of subinterval [a,b].
    h <- b - a
    c <- (a + b)/2
    d <- (a + c)/2
    e <- (c + b)/2
    fd <- f(d)
    fe <- f(e)
    fcnt <- fcnt + 2

    Q1 <- (h/6) *(fa + 4*fc + fb)               # Three point Simpson's rule
    Q2 <- (h/12)*(fa + 4*fd + 2*fc + 4*fe + fb) # Five point double Simpson's rule
    Q  <- Q2 + (Q2 - Q1)/15                     # One step of Romberg extrapolation

    if (!all(is.finite(Q)))
        stop("Improper function values: infinite or NaN encountered.")
    if (fcnt > maxfcnt)
        stop("Maximum function count exceeded; singularity likely.")
    if (abs(h) < hmin || c == a || c == b)
        stop("Minimum step size reached; singularity possible.")
    if (max(Q2 - Q) < tol)
        return(list(Q = Q, fcnt = fcnt))

    Iac <- .quadvstep(f, a, c, fa, fd, fc, tol, fcnt, hmin)
    Qac <- Iac$Q
    fcnt <- Iac$fcnt
    Icb <- .quadvstep(f, c, b, fc, fe, fb, tol, fcnt, hmin)
    Qcb <- Icb$Q
    fcnt <- Icb$fcnt

    Q <- Qac + Qcb
    return(list(Q = Q, fcnt = fcnt))
}


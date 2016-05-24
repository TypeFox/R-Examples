##
##  n e w m a r k . R  Newmark Method
##


newmark <- function(f, t0, t1, y0, ..., N = 100, zeta = 0.25, theta = 0.5) {
    stopifnot(is.numeric(y0), is.numeric(t0), length(t0) == 1, 
              is.numeric(N),  is.numeric(t1), length(t1) == 1)
    if (length(y0) != 2)
        stop("Argument 'y0' must be a numeric vector of length 2.")
    N <- floor(N)
    if (N < 2) stop("Argument 'N' must be an integer greater than 1.")

    fun <- match.fun(f)
    f <- function(t, y) fun(t, y, ...)
    #if (length(f(t0, y0)) != 2)
    #    stop("Function 'f' must always return a vector of length 2.")

    yout <- matrix(NA, N, 2)
    yout[1, ] <- y <- c(y0)

    h <- (t1 - t0)/(N-1)
    t <- 0
    ts <- linspace(t0, t1, N)

    f1 <- f(ts[1], y)

    # internal function used for root finding
    nmfun <- function(w) {
        f2 <- f(t, w)
        z1 <- w[1] - y[1] - h * y[2] - h^2 * (zeta*f2 + (0.5 - zeta) * f1)
        z2 <- w[2] - y[2] - h * ((1-theta)*f1 + theta*f2)
        c(z1, z2)
    }

    for (i in 2:N) {
        t <- ts[i]
        w <- fsolve(nmfun, y)$x
        f1 <- f(t, w)
        yout[i, ] <- w
        y <- w
    }

    return(list(t = ts, y = yout))
}

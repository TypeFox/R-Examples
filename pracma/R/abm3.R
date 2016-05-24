abm3pc <- function(f, a, b, y0, n = 50, ...) {
    stopifnot(is.numeric(a), is.numeric(b))
    stopifnot(is.numeric(y0), length(y0) == 1)
    if (!is.numeric(n) || length(n) != 1 || n < 5)
        stop("Argument 'n' must be an integer greater or equal to 5.")

    n <- floor(n)
    fun <- match.fun(f)
    f <- function(x, y) fun(x, y, ...)

    h <- (b-a)/n	
    k <- h/12

    x  <- seq(a, b, by = h)
    z  <- y <- numeric(n+1)
    z[1] <- f(a, y0)
    y[1] <- y0

    # Use midpoint method to start
    k1 <- h * z[1]
    k2 <- h * f(a + h/2, y0 + k1/2)
    k3 <- h * f(a + 0.75*h, y0 + 0.75*k2)
    y[2] <- y0 + (2*k1 +3*k2 + 4*k3)/9
    z[2] <- f(x[2], y[2])

    # Use Runge-Kutta for next step
    k1 <- h * z[2] 
    k2 <- h * f(x[2] + h/2, y[2] + k1/2)
    k3 <- h * f(x[2] + 0.75*h, y[2] + 0.75*k2)
    y[3] <- y[2] + (2*k1 +3*k2 + 4*k3)/9
    z[3] <- f(x[2], y[2])

    zz <- yy <- numeric(n)
    errorest <- numeric(n)

    # Use 3rd order A-B-M method for the remaining points
    # yy is the predicted, y the corrected value
    for (i in 3:n) {
       yy[i+1] <- y[i] + k * (23*z[i] - 16*z[i-1] + 5*z[i-2])
       zz[i+1] <- f(x[i+1], yy[i+1])
       y[i+1]  <- y[i] + k * (5*zz[i+1] + 8*z[i] - z[i-1])
       z[i+1]  <- f(x[i+1], y[i+1])

       # Error estimation
       errorest[i+1] <- -0.1 * (y[i+1] - yy[i+1])
    }

    errorest <- sqrt(abs(errorest))
    return(list(x = x, y = y, est.error = errorest))
}

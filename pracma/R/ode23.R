##
##  o d e 2 3 . R  ODE Solver
##


ode23 <- function(f, t0, tfinal, y0, ..., rtol = 1e-3, atol = 1e-6) {
    stopifnot(is.numeric(y0), is.numeric(t0), length(t0) == 1,
              is.numeric(tfinal), length(tfinal) == 1)
              
    if (is.vector(y0)) {
        y0 <- as.matrix(y0)
    } else if (is.matrix(y0)) {
        if (ncol(y0) != 1)
            stop("Argument 'y0' must be a vector or single column matrix.")
    }

    fun <- match.fun(f)
    f <- function(t, y) fun(t, y, ...)
    if (length(f(t0, y0)) != length(y0))
        stop("Argument 'f' does not describe a system of equations.")

    # Set initial parameters
    eps <- .Machine$double.eps  # Matlab parameters
    realmin <- 1e-100

    tdir <- sign(tfinal - t0)
    threshold <- atol / rtol
    hmax <- abs(0.1 * (tfinal-t0))

    t <- t0; tout <- t
    y <- y0; yout <- t(y)

    # Compute initial step size
    s1 <- f(t, y)
    r <- max(abs(s1 / max(abs(y), threshold))) + realmin
    h <- tdir * 0.8 * rtol^(1/3) / r

    # Main loop
    while (t != tfinal) {
        hmin <- 16 * eps * abs(t)
        if (abs(h) > hmax) {
            h <- tdir * hmax
        } else if (abs(h) < hmin) {
            h <- tdir * hmin
        }

        # Stretch the step if t is close to tfinal
        if (1.1 * abs(h) >= abs(tfinal - t))
            h <- tfinal - t

        # Attempt a step
        s2   <- f(t + h/2, y + h/2 * s1)
        s3   <- f(t + 3*h/4, y + 3*h/4 * s2)
        tnew <- t + h
        ynew <- y + h * (2*s1 + 3*s2 + 4*s3) / 9
        s4   <- f(tnew, ynew)

        # Estimate the error
        e   <-  h * (-5*s1 + 6*s2 + 8*s3 - 9*s4) / 72
        err <- max(abs(e / max(max(abs(y), abs(ynew)), threshold))) + realmin

        # Accept the solution if the estimated error is less than the tolerance
        if (err <= rtol) {
            t <- tnew
            y <- ynew
            tout <- c(tout, t)
            yout <- rbind(yout, t(y))
            s1 <- s4     # Reuse final function value to start new step.
        }

        # Compute a new step size
        h <- h * min(5, 0.8 * (rtol/err)^(1/3))

        # Exit early if step size is too small
        if (abs(h) <= hmin) {
            warning("Step size too small.")
            t <- tfinal
        }
    } # end while

    # Return results
    return(list(t = as.matrix(tout), y = yout))
}


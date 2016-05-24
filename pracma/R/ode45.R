##
##  o d e 5 4 . R  ODE Solver
##


ode45 = function(f, t0, tfinal, y0, ..., atol = 1e-6, hmax = 0.0)
{
    stopifnot(is.numeric(y0), is.numeric(t0), length(t0) == 1, 
              is.numeric(tfinal), length(tfinal) == 1)
    if (is.vector(y0)) {
        y0 <- as.matrix(y0)
    }
    else if (is.matrix(y0) && ncol(y0) != 1) {
        stop("Argument 'y0' must be a vector or single column matrix.")
    }
    fun <- match.fun(f)
    f <- function(t, y) fun(t, y, ...)

    pow <- 1/6                      # see p.91 in Ascher & Petzold
    nsteps <- 1e3 * (tfinal - t0)   # estimated number of steps
    if (hmax == 0.0)
        hmax <- (tfinal - t0)/2.5   # max stepsize

    # Define the Dormand-Prince 4(5) coefficients
    dp = matrix(c(0, 0, 0, 0, 0, 0,
                  1/5, 0, 0, 0, 0, 0,
                  3/40, 9/40, 0, 0, 0, 0,
                  44/45, -56/15, 32/9, 0, 0, 0,
                  19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,
                  9017/3168,  -355/33, 46732/5247, 49/176, -5103/18656, 0,
                  35/384, 0, 500/1113, 125/192, -2187/6784, 11/84),
                nrow = 7, ncol = 6, byrow=TRUE)
    b4 = c(5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40)
    b5 = c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0)
    cc = rowSums(dp)

    # Initialization
    t <- t0
    hmin <- (tfinal - t)/1e20
    h <- (tfinal - t)/100           # initial step size guess
    x <- as.matrix(y0)              # ensure x is a column vector

    nstates <- size(x,1)
    tout <- zeros(nsteps, 1)        # preallocating memory
    xout <- zeros(nsteps, nstates)

    nsteps_rej <- 0
    nsteps_acc <- 1
    tout[nsteps_acc] <- t           # first output time
    xout[nsteps_acc,] <- x          # first output solution = IC's

    # Main loop using Dormand-Prince 4(5) pair
    kk <- x %*% zeros(1, 7)
    kk[, 1] <- f(t, x)
    while (t < tfinal && h >= hmin) {
        if (t + h > tfinal) h <- tfinal - t
        for (j in 1:6)
            kk[, j+1] <- f(t + cc[j+1]*h, x + h*kk[, 1:j] %*% as.matrix(dp[j+1, 1:j]))

        # compute the 4th and 5th order estimates
        x4 <- x + h * (kk %*% b4)
        x5 <- x + h * (kk %*% b5)

        # estimate the errors
        gamma1 <- x5 - x4                   # local truncation error
        delta = Norm(gamma1, Inf);          # actual error
        tau = atol * max(Norm(x, Inf), 1.0)  # allowable error

        # update solution only if the error is acceptable
        if (delta < tau) {
            t = t + h
            x = x5              # "local extrapolation"

            nsteps_acc <- nsteps_acc + 1
            tout[nsteps_acc] <- t
            xout[nsteps_acc, ] <- x
            kk[, 1] <- kk[, 7]

        } else {                # unacceptable integration step
            nsteps_rej <- nsteps_rej + 1
        }

        # update the step size
        if (delta == 0.0) delta <- 1e-16
        h <- min(hmax, 0.8 * h * (tau/delta)^pow)
    }

    # Trim output
    if (nsteps_acc < nsteps) {
        tout <- tout[1:nsteps_acc]
        xout <- xout[1:nsteps_acc, ]
    }
    if (t < tfinal)
        warning("Step size grew too small.")

    return(list(t = tout, y = xout))
}


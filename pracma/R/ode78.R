##
##  o d e 7 8 . R  ODE Solver
##


ode78 = function(f, t0, tfinal, y0, ..., atol = 1e-6, hmax = 0.0)
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

    pow <- 1/8                      # see p.91 in the Ascher & Petzold
    if (hmax == 0.0)
        hmax <- (tfinal - t0)/2.5   # max stepsize

    # Define the Fehlberg (7,8) coefficients
    alpha <- as.matrix(c(2/27, 1/9, 1/6, 5/12, 0.5, 5/6, 1/6, 2/3, 1/3, 1, 0, 1))
    beta <- matrix( c(2/27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1/36, 1/12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          1/24, 0, 1/8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          5/12, 0, -25/16, 25/16, 0, 0, 0, 0, 0, 0, 0, 0, 0,
          0.05, 0, 0, 0.25, 0.2, 0, 0, 0, 0, 0, 0, 0, 0,
          -25/108, 0, 0, 125/108, -65/27, 125/54, 0, 0, 0, 0, 0, 0, 0,
          31/300, 0, 0, 0, 61/225, -2/9, 13/900, 0, 0, 0, 0, 0, 0,
          2, 0, 0, -53/6, 704/45, -107/9, 67/90, 3, 0, 0, 0, 0, 0,
          -91/108, 0, 0, 23/108, -976/135, 311/54, -19/60, 17/6, -1/12, 0, 0, 0, 0,
          2383/4100, 0, 0, -341/164, 4496/1025, -301/82, 2133/4100, 45/82, 45/164, 18/41, 0, 0, 0,
          3/205, 0, 0, 0, 0, -6/41, -3/205, -3/41, 3/41, 6/41, 0, 0, 0,
          -1777/4100, 0, 0, -341/164, 4496/1025, -289/82, 2193/4100, 51/82, 33/164, 12/41, 0, 1, 0),
        nrow = 13, ncol = 12, byrow = FALSE)
    chi <- as.matrix(c(0, 0, 0, 0, 0, 34/105, 9/35, 9/35, 9/280, 9/280, 0, 41/840, 41/840))
    psi <- as.matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, -1))

    # Initialization
    t <- t0
    hmin <- (tfinal - t)/1e20
    h <- (tfinal - t)/50            # initial step size guess
    x <- as.matrix(y0)              # ensure x is a column vector

    ff <- x %*% zeros(1, 13)
    tout <- t
    xout <- t(x)
    tau <- atol * max(Norm(x, Inf), 1)

    # Main loop using Fehlber (7,8) pair
    while (t < tfinal && h >= hmin) {
        if (t + h > tfinal) h <- tfinal - t

        ff[, 1] <- f(t, x)
        for (j in 1:12)
            ff[, j+1] <- f(t + alpha[j]*h, x + h*ff %*% as.matrix(beta[, j]))

        # estimate the error term
        gamma1 <- h * 41/840 * ff %*% psi   # local truncation error

        # estimate the error and the acceptable error
        delta <- Norm(gamma1, Inf);           # actual error
        tau <- atol * max(Norm(x, Inf), 1.0)  # allowable error

        # update solution only if the error is acceptable
        if (delta <= tau) {
            t <- t + h
            x <- x + h * ff %*% chi          # "local extrapolation"
            tout <- c(tout, t)
            xout <- rbind(xout, t(x))
        }

        # update the step size
        if (delta == 0.0) delta <- 1e-16
        h <- min(hmax, 0.8 * h * (tau/delta)^pow)
    }

    if (t < tfinal)
        warning("Step size grew too small: singularity likely.")

    return(list(t = tout, y = xout))
}


ode23s <- function(f, t0, tfinal, y0, jac = NULL, ...,
            rtol = 1e-03, atol = 1e-06, hmax = 0.0) {
    stopifnot(is.numeric(y0), is.numeric(t0), length(t0) == 1,
              is.numeric(tfinal), length(tfinal) == 1)
    if (t0 >= tfinal)
        warning("'t0 >= tfinal' may lead to incorrect behavior or results.")

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

    n <- length(y0); m <- length(f(t, y0))
    # use finite difference Jacobian
    if (is.null(jac)) {
        jac <- function(t, x) {
            jacob <- matrix(NA, m, n)
            hh <- numeric(n); heps <- 5e-06
            for (i in 1:n) {
                hh[i] <- heps
                jacob[, i] <- (f(t, x+hh) - f(t, x-hh)) / (2*heps)
                hh[i] <- 0
            }
            jacob
        }
    }

    # Set initial parameters
    d <- 1/(2 + sqrt(2))
    cc <- 1/2
    e32 <-  6 + sqrt(2)

    t <- t0
    tdir <- sign(tfinal - t)

    h <- tdir * 0.01 * (tfinal - t)
    if (hmax == 0.0)
        hmax <- 0.1 * abs(tfinal - t)
    hmin <- min(16 * eps(tfinal - t), h)

    y <- as.matrix(y0)
    tout <- c(t)
    yout <- c(y)

    # Main loop
    while (abs(t) < abs(tfinal) && hmin < abs(h)) {
        if (abs(t - tfinal) < abs(h))
            h <- tfinal - t

        J <- jac(t, y)

        # approximate time-derivative of f
        T <- (f(t + 0.01*h, y) - f(t, y)) / (0.01*h)

        # Wolfbrandt coefficient
        W <- eye(length(y0)) - h * d * J

        # modified Rosenbrock formula
        F1 <- f(t, y)
        k1 <- qr.solve(W, F1 + h * d * T)
        F2 <- f(t + cc * h, y + cc * h * k1)
        k2 <- qr.solve(W, F2 - k1) + k1

        # 2nd and 3rd order estimates
        ynew <- y + h * k2
        F3 = f(t + h, ynew)
        k3 = qr.solve(W, (F3 - e32*(k2 - F2) - 2*(k1 - F1) + h*d*T))

        # estimate error and acceptable error
        err <- h/6 * Norm(k1 - 2*k2 + k3)
        tau <- max(rtol * max(Norm(y),Norm(ynew)), atol)

        # check if new solution is acceptable
        if (err <= tau) {
            t <- t + h
            tout <- c(tout, t)
            y <- ynew
            yout <- cbind(yout, y)
            if (err == 0) err <- eps()/2
            # h <- min(hmax, 1.25*h)
        } else {
            if (h <= hmin)
                warning("ode23s: Requested step size too small!")
            # h <- max(hmin, 0.5*h)
        }
        # update the step size
        h <- tdir * min(hmax, abs(h)*0.8*(tau/err)^(1/3))
        if (abs(h) > hmax)
            h <- sign(h)*hmax
    }
    return(list(t = as.matrix(tout), y = t(yout)))
}


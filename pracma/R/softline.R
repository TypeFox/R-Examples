##
##  s o f t l i n e . R  Soft (Inexact) Line Search
##


softline <- function(x0, d0, f, g = NULL) {
    if (!is.numeric(x0) || !is.numeric(d0))
        stop("Arguments 'x0' and 'd0' must be numeric vectors.")
    if (length(x0) != length(d0))
        stop("Vectors 'x0' and 'd0' must have the same length.")

    f <- match.fun(f)
    if (!is.null(g)) {
        g <- match.fun(g)
    } else {
        g <- function(x) grad(f, x)
    }

    # STEP 1: Initialize search parameters
    tau <- 0.1;  chi <- 0.75
    rho <- 0.1;  sigma <- 0.1
    mhat <- 400; epsilon <- 1e-10

    xk <- c(x0)
    dk <- c(d0)

    m <- 0                              # no. of function calls
    f0 <- f(xk)
    gk <- g(xk); m <- m + 2
    deltaf0 <- f0

    # STEP 2: Initialize line search
    aL <- 0; aU <- 1e9                  # interval [a, b]
    fL <- f0

    dfL <- sum(gk * dk)                 # derivative at x0
    if (abs(dfL) > epsilon)   a0 <- -2*deltaf0 / dfL
    else                      a0 <- 1
    if (a0 <= 1e-9 || a0 > 1) a0 <- 1
    
    # STEP 3 and 4: Estimate a0 and compute f0
    repeat {
        deltak <- a0*dk
        f0 <- f(xk + deltak); m <- m + 1

        # STEP 5: Interpolation
        if (f0 > (fL + rho*(a0 - aL)*dfL) &&
            abs(fL - f0) > epsilon && m < mhat) {
            if (a0 < aU) aU <- a0

            # Compute a0hat by extrapolation
            a0hat <- aL + ((a0 - aL)^2 * dfL) / (2*(fL - f0 + (a0 - aL)*dfL))
            a0Lhat <- aL + tau * (aU - aL)
            if (a0hat < a0Lhat)
                a0hat <- a0Lhat
            a0Uhat <- aU - tau * (aU - aL)
            if (a0hat > a0Uhat)
                a0hat <- a0Uhat
            a0 <- a0hat

        # STEP 6: Compute df0
        } else {
            df0 <- sum(g(xk + a0*dk) * dk)
            m <- m + 1

            # STEP 7: Extrapolation
            if (df0 < sigma*dfL &&
                abs(fL - f0) > epsilon && m < mhat && dfL != df0) {
                deltaa0 <- (a0 - aL) * df0 / (dfL - df0)
                if (deltaa0 <= 0) a0hat <- 2*a0
                else              a0hat <- a0 + deltaa0

                a0Uhat <- a0 + chi * (aU - a0)
                if (a0hat > a0Uhat) a0hat <- a0Uhat

                # Pepare next iteration
                aL <- a0
                a0 <- a0hat
                fL <- f0
                dfL <- df0
            } else {
                break
            }
        }
    }
    z <- max(a0, 1e-5)
    return(z)
}

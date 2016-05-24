##
##  s p i n t e r p . R
##


spinterp <- function(x, y, xp) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(xp))
    stopifnot(is.vector(x), is.vector(y), is.vector(xp))
    n <- length(x)
    n1 <- n - 1
    if (n <= 3)
        stop("Length of arguments 'x', 'y' must be greater than 3.")
    if (length(y) != n)
        stop("Arguments 'x', 'y' must be vectors of the same length.")
    if(any(is.na(y)))
        stop("NAs are not allowed in argument 'y'.")

    # M o n o t o n i c i t y
    h <- diff(x)
    dy <- diff(y)
    if (any(h <= 0))
        stop("Argument 'x' must be a sorted list af real values.")
    if (any(dy < 0))
        stop("Argument 'y' must be monotonically increasing.")

    # C o n v e x i t y
    delta <- dy / h
    dd <- diff(delta)
    cnvx <- if (all(dd > 0)) TRUE else FALSE

    # Approximate the derivatives at all data points
    if (cnvx) d_mode <- "harmonic"
    else      d_mode <- "geometric"
    d <- rep(NA, n)
    d[1] <- delta[1]
    d[n] <- delta[n1]
    if (d_mode == "arithmetic") {
        for (j in 2:n1)
            d[j] <- (h[j]*delta[j-1] + h[j-1]*delta[j]) / (h[j] + h[j-1])
    } else if (d_mode == "geometric") {
        for (j in 2:n1) {
            d[j] <- (delta[j-1]^h[j] * delta[j]^h[j-1])^(1/(h[j-1]+h[j]))
        }
    } else if (d_mode == "harmonic") {
        for (j in 2:n1) {
            d[j] <- (h[j] + h[j-1]) / (h[j]/delta[j-1] + h[j-1]/delta[j])
        }
    }

    # "C u b i c i t y"
    r <- rep(3, n1)
    if (!cnvx) r_mode <- "monotonic"  #  a n d  monotone
    else       r_mode <- "otherwise"
    r_mode <- "monotonic"               # Fix to "monotonic" for the moment !
    # Now define the r-values for Delbourg & Gregory
    if (r_mode == "monotonic") {
        k <- which(delta != 0)
        #r <- 1 + (d[1:n1] + d[2:n]) / delta  # strictly monotonic
        r[k] <- (d[k] + d[k+1]) / delta[k]
    } else if (r_mode == "otherwise") {
        for (j in 1:n1) {
            r[j] <- 1 + (d[j+1]-delta[j]) / (delta[j]-d[j]) +
                        (delta[j]-d[j]) / (d[j+1]-delta[j])
        }
    }

    # Apply cubic Delbourgo & Gregory formula
    m <- length(xp)
    fi <- findInterval(xp, x)
    yp <- rep(NA, m)

    for (j in 1:m) {
        i <- fi[j]  # findInterval(xp[j], x)
        if (i < n) {
            theta <- (xp[j] - x[i]) / h[i]
            P <- y[i+1] * theta^3 + 
                 (r[i]*y[i+1] - h[i]*d[i+1]) * theta^2 * (1-theta) + 
                 (r[i]*y[i] + h[i]*d[i]) * theta * (1-theta)^2 + 
                 y[i] * (1-theta)^3
            Q <- 1 + (r[i]-3) * theta * (1-theta)
            yp[j] <- P / Q
        } else {
            yp[j] <- y[n]
        }
    }

    return(yp)
}

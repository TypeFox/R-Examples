##
##  g r a d i e n t . R  Discrete Derivatives
##


gradient <- function(F, h1 = 1, h2 = 1) {
    if (length(F) == 0 )
        return(c())
    if (!is.numeric(F))
        stop("Argument 'F' must be a numeric vector or matrix.")
    if (length(h1) == 0 || length(h2) == 0 ||
        (length(h1) == 1 && h1 == 0) || (length(h2) == 1 && h2 == 0))
        stop("Arguments 'h1', 'h2' must be non-empty and non-zero.")
    if (any(diff(h1) == 0) || any(diff(h2) == 0))
        stop("Arguments 'h1' and 'h2' must be strictly increasing.")

    if (is.vector(F)) {
        n <- length(F)
        if (n == 1) return(0)
        if (length(h1) == 1) {
            x <- seq(1*h1, n*h1, length.out = n)
        } else if (length(h1) == n) {
            x <- h1
        } else
            stop("Length of 'h1' must be 1 or equal to length of 'F'.")

        g <- numeric(n)
        g[1] <- (F[2] - F[1]) / (x[2] - x[1])
        g[n] <- (F[n] - F[n-1]) / (x[n] - x[n-1])

        if (n > 2)
            g[2:(n-1)] <- (F[3:n] - F[1:(n-2)]) / (x[3:n] - x[1:(n-2)])

        return(g)
        
    } else if (is.matrix(F)) {
        # stop("Two-dimensional version not yet implemented.")
        n <- nrow(F)
        m <- ncol(F)

        if (length(h1) == 1) {
            x <- seq(1*h1, m*h1, length.out = m)
        } else if (length(h1) == m) {
            x <- h1
        } else
            stop("Length of 'h1' must be 1 or equal to ncol of 'F'.")
        if (length(h2) == 1) {
            y <- seq(1*h2, n*h2, length.out = n)
        } else if (length(h2) == n) {
            y <- h2
        } else
            stop("Length of 'h2' must be 1 or equal to nrow of 'F'.")

        gX <- gY <- 0 * F  # matrix(NA, nrow = n, ncol = m)

        gX[, 1] <- (F[, 2] - F[, 1]) / (x[2] - x[1])
        gX[, m] <- (F[, m] - F[, m-1]) / (x[m] - x[m-1])
        if (m > 2)
            gX[, 2:(m-1)] <- (F[, 3:m] - F[, 1:(m-2)]) / (x[3:m] - x[1:(m-2)])

        gY[1, ] <- (F[2, ] - F[1, ]) / (y[2] - y[1])
        gY[n, ] <- (F[n, ] - F[n-1, ]) / (y[n] - y[n-1])
        if (n > 2)
            gY[2:(n-1), ] <- (F[3:n, ] - F[1:(n-2), ]) / (y[3:n] - y[1:(n-2)])

        return(list(X = gX, Y = gY))

    } else
        stop("Argument 'F' cannot be a higher-dimensional array.")
}

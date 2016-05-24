ratinterp <- function(x, y, xs = x) {
    if (!is.vector(x, mode="numeric") || !is.vector(y, mode="numeric"))
        stop("Arguments 'x' and 'y' must be numeric vectors.")
    m <- length(x)
    if (length(y) != m)
        stop("Arguments 'x' and 'y' must be vectors of the same length.")
    if (m <= 2)
        stop("Arguments 'x', 'y' must have at least a length >= 3.")

    if (is.unsorted(x))
        stop("Argument 'x' must be a sorted vector")
    
    n <- length(xs)
    ys <- numeric(n)
    node_p <- FALSE

    for (h in 1:n) {
        for (i in 1:m) {
            d <- xs[h] - x[i]
            # First stage: xs[h] is data points
            if (d == 0) {
                ys[h] <- y[i]
                node_p <- TRUE
            }
        }

        if (!node_p) {
            # Second stage
            R <- matrix(0, nrow = m, ncol = m)
            R[, 1] <- y
            for (i in 1:(m-1)) {
                D <- R[i+1, 1] - R[i, 1]
                rr <- (xs[h] - x[i]) / (xs[h] - x[i+1])
                denom <- rr * (1 - D/R[i+1, 1]) - 1
                R[i, 2] <- R[i+1, 1] + D/denom
            }

            for (j in 3:m) {
                # Third and next stages
                for (i in 1:(m-j+1)) {
                    D <- R[i+1, j-1] - R[i, j-1]
                    rr <- (xs[h] - x[i]) / (xs[h]- x[i+j-1])
                    if (D == 0) {
                        R[i, j] <- R[i+1, j-1]
                    } else {
                        DD <- R[i+1, j-1] - R[i+1, j-2]
                        denom <- rr * (DD - D) - DD
                        R[i, j] <- R[i+1, j-1] + D * DD / denom
                    }
                }
            }
            ys[h] <- R[1, m]
        }
        node_p <- FALSE
    }
    return(ys)
}

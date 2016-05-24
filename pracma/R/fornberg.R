fornberg <- function(x, y, xs, k = 1) {
    stopifnot(is.numeric(x), is.numeric(y), is.numeric(xs))
    if (any(is.na(y))) {
        inna <- which(!is.na(y))
        x <- x[inna]; y <- y[inna]
    }
    n <- length(x); l <- length(xs)
    if (length(unique(x)) != n)
        stop("All elements in vector 'x' must be different.")
    if (k >= n)
        stop("Length of 'x' must be greater than k.")
    if (k <= 0)
        stop("Order 'k' must be between 1 and length of 'x'.")

    m <- k
    Y <- matrix(NA, nrow = l, ncol = k+1)

    for (ij in 1:l) {
        x0 <- xs[ij]
        c1 <- 1
        c4 <- x[1] - x0
        C <- zeros(n, m+1)
        C[1,1] <- 1
    
        for (i in 1:(n-1)) {
            i1 <- i+1
            mn <- min(i,m)
            c2 <- 1
            c5 <- c4
            c4 <- x[i1] - x0
            for (j in 0:(i-1)) {
                j1 <- j+1
                c3 <- x[i1] - x[j1]
                c2 <- c2*c3
                if (j == i-1) {
                    for (s in mn:1) {
                        s1 <- s+1
                        C[i1,s1] <- c1*(s*C[i1-1,s1-1] - c5*C[i1-1,s1])/c2
                    }
                    C[i1,1] <- -c1*c5*C[i1-1,1]/c2
                }
                for (s in mn:1) {
                    s1 <- s+1
                    C[j1,s1] <- (c4*C[j1,s1] - s*C[j1,s1-1])/c3
                }
                C[j1,1] <- c4*C[j1,1]/c3
            }
            c1 <- c2
        }
        Y[ij, ] <- y %*% C
    }

    if (k == 0) Y <- Y[, 1, drop = FALSE]
    return(Y)
}

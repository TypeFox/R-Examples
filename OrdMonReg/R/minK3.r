minK3 <- function(f, w, n){
    x <- rep(NA, 2 * n)
    i <- 1:n

    ind1 <- (1:n)[(f[i] <= f[n + i])]
    ind2 <- (1:n)[(f[i] > f[n + i])]

    ## no violation
    x[ind1] <- f[ind1]
    x[n + ind1] <- f[n + ind1]

    ## violation
    x[ind2] <- (w[ind2] * f[ind2] + w[n + ind2] * f[n + ind2]) / (w[ind2] + w[n + ind2])
    x[n + ind2] <- x[ind2]

    return(x)
    }

## cyq.R -- Fletcher's Chebyquad problem

rm(list = ls())
library(Rcgmin)
# Ref: Fletcher, R. (1965) Function minimization without
#   calculating derivatives -- a review,
#         Computer J., 8, 33-41.

# Note we do not have all components here e.g., .jsd, .h

cyq.f <- function(x) {
    rv <- cyq.res(x)
    f <- sum(rv * rv)
}

cyq.res <- function(x) {
    # Fletcher's chebyquad function m = n -- residuals
    n <- length(x)
    res <- rep(0, n)  # initialize
    for (i in 1:n) {
        #loop over resids
        rr <- 0
        for (k in 1:n) {
            z7 <- 1
            z2 <- 2 * x[k] - 1
            z8 <- z2
            j <- 1
            while (j < i) {
                z6 <- z7
                z7 <- z8
                z8 <- 2 * z2 * z7 - z6  # recurrence to compute Chebyshev polynomial
                j <- j + 1
            }  # end recurrence loop
            rr <- rr + z8
        }  # end loop on k
        rr <- rr/n
        if (2 * trunc(i/2) == i) {
            rr <- rr + 1/(i * i - 1)
        }
        res[i] <- rr
    }  # end loop on i
    res
}

cyq.jac <- function(x) {
    #  Chebyquad Jacobian matrix
    n <- length(x)
    cj <- matrix(0, n, n)
    for (i in 1:n) {
        # loop over rows
        for (k in 1:n) {
            # loop over columns (parameters)
            z5 <- 0
            cj[i, k] <- 2
            z8 <- 2 * x[k] - 1
            z2 <- z8
            z7 <- 1
            j <- 1
            while (j < i) {
                # recurrence loop
                z4 <- z5
                z5 <- cj[i, k]
                cj[i, k] <- 4 * z8 + 2 * z2 * z5 - z4
                z6 <- z7
                z7 <- z8
                z8 <- 2 * z2 * z7 - z6
                j <- j + 1
            }  # end recurrence loop
            cj[i, k] <- cj[i, k]/n
        }  # end loop on k
    }  # end loop on i
    cj
}


cyq.g <- function(x) {
    cj <- cyq.jac(x)
    rv <- cyq.res(x)
    gg <- as.vector(2 * rv %*% cj)
}

cat("Fletcher chebyquad function in file cyq.R\n")

nn <- c(2, 3, 5, 8, 10, 20, 30)

for (n in nn) {
    cat("Chebyquad in ", n, " parameters\n")
    afname <- paste("acyq", n, "G.txt", sep = "")
    lower <- rep(-10, n)
    upper <- rep(10, n)
    bdmsk <- rep(1, n)  # free all parameters
    x0 <- 1:n
    x0 <- x0/(n + 1)  # Initial value suggested by Fletcher
    ut <- system.time(ans <- Rcgmin(x0, cyq.f, cyq.g, lower, 
        upper, bdmsk, control = list(trace = 1)))[1]
    print(ans)
    cat("time = ", ut, "\n")
}

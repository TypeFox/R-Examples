## sim.coalescent.R (2014-07-15)

##   Coalescent Simulation and Visualisation

## Copyright 2014 Emmanuel Paradis

## This file is part of the R-package `coalescentMCMC'.
## See the file ../COPYING for licensing issues.

sim.coalescent <-
    function(n = 5, TIME = 50, growth.rate = NULL, N.0 = 50, N.final = 20,
             col.lin = "grey", col.coal = "blue", pch = NULL, ...)
{
    if (n > N.final) stop("sample size smaller than final population size")
    if (is.null(TIME)) TIME <- log(N.final/N.0) / growth.rate
    if (is.null(growth.rate)) growth.rate <- log(N.final/N.0) / TIME
    if (is.null(N.0)) N.0 <- N.final/exp(growth.rate * TIME)
    if (is.null(N.final)) N.final <- N.0 * exp(growth.rate * TIME)

    N <- round(N.0 * exp(growth.rate * 1:TIME))
    x <- numeric(0)
    for (j in 1:TIME) x <- c(x, (-N[j]/2 + 0.5):(N[j]/2 - 0.5))
    y <- rep(1:TIME, N)
    x1 <- x2 <- y1 <- y2 <- numeric(0)

    for (j in TIME:2) {
        anc <- sort(sample(1:N[j - 1], N[j], replace = TRUE))
        x1 <- c(x1, x[which(y == j)])
        y1 <- c(y1, rep(j, N[j]))
        x2 <- c(x2, x[which(y == j - 1)][anc])
        y2 <- c(y2, rep(j - 1, N[j]))
    }

    plot(x, y, type = "n", xlab = "", ylab = "", frame = FALSE,
         xaxt = "n", yaxt = "n")
    segments(x1, y1, x2, y2, col = col.lin)

    dat <- sample(x[which(y == TIME)], n)
    for (j in TIME:2) {
        for (i in seq_len(n)) {
            index <- which(x1 == dat[i] & y1 == j)
            segments(x1[index], j, x2[index], j - 1, lwd = 2, col = col.coal)
            dat[i] <- x2[index]
        }
    }
    if (!is.null(pch))
        points(x, y, pch = pch, ...)
}

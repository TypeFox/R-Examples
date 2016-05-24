## R2.test.R (2009-05-13)

##   Ramos-Onsins--Rozas Test of Neutrality

## Copyright 2009 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

R2.test <- function(x, B = 1000, theta = 1, plot = TRUE, quiet = FALSE, ...)
{
    if (is.list(x)) x <- as.matrix(x)
    n <- dim(x)[1]
    k <- mean(dist.dna(x, "N"))
    ss <- seg.sites(x)
    U <- numeric(n)
    for (i in 1:n) for (j in ss)
        if (all(x[i, j, drop = TRUE] != x[-i, j])) U[i] <- U[i] + 1
    U <- (U - k/2)^2
    R2.obs <- sqrt(sum(U)/n)/length(ss)
    if (B) {
        R <- numeric(B)
        if (!quiet) progbar <- utils::txtProgressBar(style = 3)
        for (b in 1:B) {
            if (!quiet) utils::setTxtProgressBar(progbar, b/B)
            tr <- rcoal(n, rep("", n))
            tr$edge.length <- rpois(2*n - 2, theta * tr$edge.length)
            d <- cophenetic(tr)
            k <- mean(d)
            U <- tr$edge.length[tr$edge[, 2] <= n]
            U <- (U - k/2)^2
            R[b] <- sqrt(sum(U)/n)/sum(tr$edge.length)
        }
        if (!quiet) close(progbar)
        if (plot) {
            hist(R, xlab = expression("Simulated "*italic(R)[2]),
                 main = "", ...)
            abline(v = R2.obs, lty = 2)
            mtext(expression("Calculated "*italic(R)[2]), at = R2.obs)
        }
        R <- na.omit(R)
        list(R2 = R2.obs, P.val = sum(R < R2.obs)/length(R))
    } else R2.obs
}

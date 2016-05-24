# exampleLS2.R -- version 2011-01-11
# ... continues exampleLS.R
trials <- 1000L

lambda <- 1; betaTRUE <- c(4,-2,2,lambda)
allResults <- array(NA, dim = c(trials, 3L))
for (t in 1L:trials) {
    yM <- NS(betaTRUE,tm) + rnorm(length(tm), mean = 0, sd = 0.01)
    allResults[t, ] <- qr.solve(NSf(lambda, tm), yM)
}

# compare results
par(mfrow = c(1L, 3L))
for (i in 1L:3L) plot(ecdf(allResults[ ,i]))

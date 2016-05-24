# exampleLS.R -- version 2010-12-18
# set up artificial yield curve and plot it
tm <- 1L:10L; paramTRUE <- c(4,-2,2,1)
yM <- NS(paramTRUE,tm)
plot(tm, yM, xlab = "maturities in years",
             ylab = "yields in %-points")

# fix lambda and run regression
lambda <- 1.5
result <- lm(yM ~ -1 + NSf(lambda,tm))

# compare results
plot(yM - result$fitted.values, xlab = "maturities in years",
                                ylab = "errors in %-points")

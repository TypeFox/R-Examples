## -*- truncate-lines: t; -*-

test.ytm <- function() {
    cf <- c(5, 5, 5, 5, 5, 105)   ## cashflows
    times <- 1:6                  ## maturities
    y <- 0.0127                   ## the "true" yield
    b0 <- vanillaBond(cf, times, yields = y)
    cf <- c(-b0, cf); times <- c(0, times)    
    checkEquals(y, ytm(cf, times), tolerance = 1e-4)
    checkException(checkEquals(y, ytm(cf, times), tolerance = 1e-7),
                   silent = TRUE)
    checkEquals(y, ytm(cf, times, tol = 1e-8), tolerance = 1e-7)


    cf <- c(5, 5, 5, 5, 5, 105)   ## cashflows
    times <- 1:6                  ## maturities
    y <- 0.02+0.01                ## the "true" yield: 2% risk-free+1% premium
    b0 <- vanillaBond(cf, times, yields = y)
    cf <- c(-b0, cf); times <- c(0, times)    
    checkEquals(ytm(cf, times)-ytm(cf, times, offset = 0.02), 0.02, tolerance=1e-5)

    cf <- c(5, 5, 5, 5, 5, 105)   ## cashflows
    times <- 1:6                  ## maturities
    y <- NS(c(6,9,10,5)/100, times) ## 1%premium
    b0 <- vanillaBond(cf, times, yields = y + 0.01)
    cf <- c(-b0, cf); times <- c(0, times)    
    checkTrue(abs(ytm(cf, times, offset = c(0,y)) - 0.01) < 1e-5)
}


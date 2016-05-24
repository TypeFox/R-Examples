library(forecast)
library(tseries)

test_that('frbe', {
    res <- frbe(as.ts(1:10 + runif(10)), h=4)
    
    expect_true(is.frbe(res))
})

context("Fast.approx")

test_that("fast reshape I", {
    set.seed(1)
    x <- sort(rnorm(1e5))
    y <- rnorm(1e3)
    val <- x[fast.approx(x,y)]
    val2 <- sapply(y, function(y) x[which.min(abs(x-y))])
    expect_identical(val,val2)
    ##rbenchmark::benchmark(fast.approx(x,y,type="left"),replications=1000)
    ##rbenchmark::benchmark(prodlim::sindex(x,y),replications=1000)    
    val <- fast.approx(x,y,type="left") # Number of observations in x less than y
    val3 <- prodlim::sindex(x,y)
    expect_identical(val,val3)    
})

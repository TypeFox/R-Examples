context("Kernel")

test_that("Exponential kernel", {
    X <- matrix(1:4, ncol=4)
    Xa <- matrix(1:4+1, ncol=4)
    sigma <- 2
    expect_that( pattern(X, Xa, sigma), equals(exp(-0.5)) )
})

test_that("Apply kernel over all patterns", {
    X <- matrix(1:4, ncol=4)
    Xa <- matrix(rep(1:4+1,5), ncol=4, byrow=TRUE)
    sigma <- 2
    expect_that( patterns(Xa, X, sigma), equals(rep(exp(-0.5),5)) )
})

test_that("Compute the probability density function", {
    X <- matrix(1:4, ncol=4)
    Xa <- matrix(rep(1:4+1,5), ncol=4, byrow=TRUE)
    sigma <- 2
    expected <- 1 /((2 * pi) ^ (4 / 2) * sigma ^ 4) / 5 * exp(-0.5) * 5
    expect_that( fA(Xa, X, sigma), equals(expected) )
})
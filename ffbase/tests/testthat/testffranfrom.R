
library(testthat)
library(ff)


test_that("ffrandom works",{
    set.seed(123)
    rnorm(1)
    a <- rnorm(100)
    set.seed(123)
    b <- ffrandom(100, rnorm)
    expect_equal(a, b[])
})

library("testthat")
library("permute")

context("Testing nobs() methods")

test_that("numeric nobs method", {
    v <- runif(10)
    n <- nobs(v)
    expect_identical(n, 10L)
})

test_that("integer nobs method", {
    v <- 1L:10L
    n <- nobs(v)
    expect_identical(n, 10L)
})

test_that("matrix nobs method", {
    m <- matrix(1:9, nrow = 3)
    n <- nobs(m)
    expect_identical(n, 3L)
})

test_that("data frame nobs method", {
    df <- as.data.frame(matrix(1:9, nrow = 3))
    n <- nobs(df)
    expect_identical(n, 3L)
})

test_that("factor nobs method", {
    f <- factor(c(1,2,3,2,1,4,5,2,1,4))
    n <- nobs(f)
    expect_identical(n, 10L)
})

library("testthat")
library("permute")

context("Testing blocks() and Blocks()")

test_that("Blocks() returns a list with 'strata'", {
    f <- gl(4, 10)
    b <- Blocks(f)
    expect_is(b, "list")
    expect_true("strata" %in% names(b))
})

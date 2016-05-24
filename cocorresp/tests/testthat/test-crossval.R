## Test crossval() method

library("testthat")
library("cocorresp")

context("Testing crossval() method")

## load data
data(beetles, plants, package = "cocorresp")
beetles <- log(beetles + 1)            # log transform the bettle data

test_that("crossval() works", {
    skip_on_cran()
    bp.loo <- crossval(beetles, plants)
    expect_is(bp.loo, "crossval")
    expect_is(bp.loo, "list")
    expect_named(bp.loo, c("dimx","dimy","n.axes","press0",
                           "CVfit","varianceExp","totalVar", "call",
                           "nam.dat"))
    expect_output(print(bp.loo),
                  regexp = "Cross-validation for Predictive Co-Correspondence Analysis")
})

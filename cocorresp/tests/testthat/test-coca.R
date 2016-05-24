## Test `coca()` function

## load packages
library("testthat")
library("cocorresp")

context("Testing coca()")

## Load data
data(beetles)
beetles <- log(beetles + 1)            # log transform the bettle data
data(plants)
data(bryophyte)
data(vascular)

test_that("coca() works & returns correct object", {
    ## symmetric CoCA
    expect_message(bp.sym <- coca(beetles ~ ., data = plants, method = "symmetric"),
                   regexp = "some species contain no data")
    expect_is(bp.sym, "coca")
    expect_is(bp.sym, "symcoca")
    expect_is(bp.sym, "list")
    expect_named(bp.sym, c("scores", "lambda", "X", "loadings",
                           "residuals", "inertia", "rowsum", "colsum",
                           "nam.dat", "n.axes", "weights", "call"))
    expect_is(bp.sym$scores, "list")
    expect_output(print(bp.sym), regexp = "Symmetric Co-Correspondence Analysis")
    expect_silent(plot(bp.sym))

    ## predictive CoCA
    expect_silent(m <- coca(y = bryophyte, x = vascular,
                             reg.method = "eigen"))
    expect_is(m, "coca")
    expect_is(m, "predcoca")
    expect_is(m, "list")
    expect_named(m, c("nam.dat", "method", "n.axes", "lambda",
                      "scores", "R0", "Ychi", "call",
                      "loadings", "fitted", "varianceExp",
                      "totalVar"))
    expect_is(m$scores, "list")
    expect_output(print(m), regexp = "Predictive Co-Correspondence Analysis")
    expect_silent(plot(m))

    ## formula and SIMPLS
    expect_silent(m2 <- coca(bryophyte ~ ., data = vascular, reg.method = "simpls"))
    expect_is(m2, "coca")
    expect_is(m2, "predcoca")
    expect_is(m2, "list")
    expect_named(m2, c("nam.dat", "call", "method", "scores", "loadings", "fitted",
                       "varianceExp", "totalVar", "lambda", "n.axes", "Ychi",
                       "R0"))
    expect_is(m2$scores, "list")
    expect_output(print(m2), regexp = "Predictive Co-Correspondence Analysis")
    expect_silent(plot(m2))
})

## load packages
library("testthat")
library("cocorresp")

context("Testing envfit()")

## symmetric CoCA
data(beetles, plants, verges)
beetles <- log(beetles + 1)            # log transform the bettle data


test_that("envfit() method works for symmetric CoCA", {
    ## fit the model
    expect_message(bp.sym <- coca(beetles ~ ., data = plants,
                                  method = "symmetric"),
                   regexp = "some species")
    ## fit vectors for the environmental data
    sol <- envfit(bp.sym, verges, which = "response")
    expect_output(print(sol), regexp = "VECTORS")
    expect_is(sol, "envfit")

    ## plot the response matrix and the fitted vectors
    expect_silent(plot(bp.sym, which = "response"))
    expect_silent(plot(sol, add = TRUE))
})

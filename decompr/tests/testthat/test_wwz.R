# load the package
library(decompr)

# load test data
data(leather)

# leontief decomposition
w <- decomp(inter,
            final,
            countries,
            industries,
            out,
            method = "wwz")

# define context
context("output format")

test_that("output size matches", {
  expect_equal(length(w), 29 )
  expect_equal(dim(w)[1], 27 )
})

test_that("output format matches", {
  expect_match(typeof(w[,5]), "double")
})
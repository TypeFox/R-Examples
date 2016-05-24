stopifnot(require(testthat))
context("issue#2")

library(naturalsort)

test_that("Produces an error when all elements are empty character", {
   text <- ""
   expected <- ""
   actual <- naturalsort(text)
   expect_that(actual, equals(expected))
})

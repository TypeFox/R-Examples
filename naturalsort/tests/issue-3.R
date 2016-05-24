stopifnot(require(testthat))
context("issue#3")

library(naturalsort)

test_that("Produces an error when input is empty", {
   text <- character(0)
   expected <- character(0)
   actual <- naturalsort(text)
   expect_that(actual, equals(expected))
})

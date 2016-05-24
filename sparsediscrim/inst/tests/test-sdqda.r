library(testthat)
library(sparsediscrim)

context("The SDQDA Classifier from Pang et al. (2009)")

test_that("The SDQDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  sdqda_out <- sdqda(Species ~ ., data = iris[train, ])
  predicted <- predict(sdqda_out, iris[-train, -5])$class

  sdqda_out2 <- sdqda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(sdqda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the SDQDA classifier
  expect_equal(predicted, predicted2)
})

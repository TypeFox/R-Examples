library(testthat)
library(sparsediscrim)

context("The SDLDA Classifier from Pang et al. (2009)")

test_that("The SDLDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  sdlda_out <- sdlda(Species ~ ., data = iris[train, ])
  predicted <- predict(sdlda_out, iris[-train, -5])$class

  sdlda_out2 <- sdlda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(sdlda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the SDLDA classifier
  expect_equal(predicted, predicted2)
})

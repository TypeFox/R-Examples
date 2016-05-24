library(testthat)
library(sparsediscrim)

context("The DQDA Classifier from Dudoit et al. (2002)")

test_that("The DQDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  dqda_out <- dqda(Species ~ ., data = iris[train, ])
  predicted <- predict(dqda_out, iris[-train, -5])$class

  dqda_out2 <- dqda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(dqda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the DQDA classifier
  expect_equal(predicted, predicted2)
})

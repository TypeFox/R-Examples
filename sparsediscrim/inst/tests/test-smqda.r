library(testthat)
library(sparsediscrim)

context("The SmDQDA Classifier from Tong et al. (2012)")

test_that("The SmDQDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  smdqda_out <- smdqda(Species ~ ., data = iris[train, ])
  predicted <- predict(smdqda_out, iris[-train, -5])$class

  smdqda_out2 <- smdqda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(smdqda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the SmDQDA classifier
  expect_equal(predicted, predicted2)
})

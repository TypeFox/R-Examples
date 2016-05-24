library(testthat)
library(sparsediscrim)

context("The SmDLDA Classifier from Tong et al. (2012)")

test_that("The SmDLDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  smdlda_out <- smdlda(Species ~ ., data = iris[train, ])
  predicted <- predict(smdlda_out, iris[-train, -5])$class

  smdlda_out2 <- smdlda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(smdlda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the SmDLDA classifier
  expect_equal(predicted, predicted2)
})

library(testthat)
library(sparsediscrim)

context("The DLDA Classifier from Dudoit et al. (2002)")

test_that("The DLDA classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  dlda_out <- dlda(Species ~ ., data = iris[train, ])
  predicted <- predict(dlda_out, iris[-train, -5])$class

  dlda_out2 <- dlda(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(dlda_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the DLDA classifier
  expect_equal(predicted, predicted2)
})

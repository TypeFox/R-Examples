library(testthat)
library(sparsediscrim)

context("The MDMP Classifier from Srivistava and Kubokawa (2007)")

test_that("The MDMP classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  mdmp_out <- mdmp(Species ~ ., data = iris[train, ])
  predicted <- predict(mdmp_out, iris[-train, -5])$class

  mdmp_out2 <- mdmp(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(mdmp_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the MDMP classifier
  expect_equal(predicted, predicted2)
})

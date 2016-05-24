library(testthat)
library(sparsediscrim)

context("The MDEB Classifier from Srivistava and Kubokawa (2007)")

test_that("The MDEB classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  mdeb_out <- mdeb(Species ~ ., data = iris[train, ])
  predicted <- predict(mdeb_out, iris[-train, -5])$class

  mdeb_out2 <- mdeb(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(mdeb_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the MDEB classifier
  expect_equal(predicted, predicted2)
})

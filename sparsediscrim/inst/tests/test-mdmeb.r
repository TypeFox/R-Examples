library(testthat)
library(sparsediscrim)

context("The MDMEB Classifier from Srivistava and Kubokawa (2007)")

test_that("The MDMEB classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  mdmeb_out <- mdmeb(Species ~ ., data = iris[train, ])
  predicted <- predict(mdmeb_out, iris[-train, -5])$class

  mdmeb_out2 <- mdmeb(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(mdmeb_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the MDMEB classifier
  expect_equal(predicted, predicted2)
})

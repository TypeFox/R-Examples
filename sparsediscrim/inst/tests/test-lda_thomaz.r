library(testthat)
library(sparsediscrim)

context("LDA with the Thomaz-Kitani-Gillies Covariance Matrix")

test_that("The lda_thomaz classifier works properly on the iris data set", {
  require('MASS')

  set.seed(42)
  n <- nrow(iris)
  train <- sample(seq_len(n), n / 2)
  lda_thomaz_out <- lda_thomaz(Species ~ ., data = iris[train, ])
  predicted <- predict(lda_thomaz_out, iris[-train, -5])$class

  lda_thomaz_out2 <- lda_thomaz(x = iris[train, -5], y = iris[train, 5])
  predicted2 <- predict(lda_thomaz_out2, iris[-train, -5])$class

  # Tests that the same labels result from the matrix and formula versions of
  # the lda_thomaz classifier
  expect_equal(predicted, predicted2)
})

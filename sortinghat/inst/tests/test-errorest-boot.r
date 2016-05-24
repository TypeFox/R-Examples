library('testthat')
library('sortinghat')

context("Bootstrap Error Rate")

test_that("Error rate works correctly on the Iris data set using MASS:::lda", {
  require('MASS')
  iris_x <- data.matrix(iris[, -5])
  iris_y <- iris[, 5]

  lda_wrapper <- function(object, newdata) {
    predict(object, newdata)$class
  }

  set.seed(42)
  error_rate <- errorest_boot(x = iris_x, y = iris_y, train = MASS:::lda,
                              classify = lda_wrapper)

  # This value was computed previously.
  expected_estimate <- 0.0228

  expect_equal(error_rate, expected_estimate, tolerance = 1e-6)
})


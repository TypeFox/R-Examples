library('testthat')
library('sortinghat')

context("Bootstrap Cross-Validation Error Rate")

test_that("Error rate works correctly on the Iris data set using MASS:::lda", {
  require('MASS')
  iris_x <- data.matrix(iris[, -5])
  iris_y <- iris[, 5]

  lda_wrapper <- function(object, newdata) {
    predict(object, newdata)$class
  }

  # Testing default arguments
  set.seed(42)
  error_rate <- errorest_bcv(x = iris_x, y = iris_y, train = MASS:::lda,
                             classify = lda_wrapper)

  # This value was computed previously.
  expected_estimate <- 0.02213333  

  expect_equal(error_rate, expected_estimate, tolerance = 1e-6)

  # Testing some arguments
  set.seed(42)
  error_rate <- errorest_bcv(x = iris_x, y = iris_y, train = MASS:::lda,
                             classify = lda_wrapper, num_bootstraps = 10,
                             hold_out = 1)

  # This value was computed previously.
  expected_estimate <- 0.01866667

  expect_equal(error_rate, expected_estimate, tolerance = 1e-6)
})


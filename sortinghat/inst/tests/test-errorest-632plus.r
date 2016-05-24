library('testthat')
library('sortinghat')

context(".632+ Error Rate")

test_that("Error rate works correctly on the Iris data set using MASS:::lda", {
  require('MASS')
  iris_x <- data.matrix(iris[, -5])
  iris_y <- iris[, 5]

  lda_wrapper <- function(object, newdata) {
    predict(object, newdata)$class
  }

  # Computes the apparent and LOO-Boot error rates up front b/c there are
  # instances where these are precomputed (e.g., simulations)
  set.seed(42)
  apparent <- errorest_apparent(x = iris_x, y = iris_y, train = MASS:::lda,
                                classify = lda_wrapper)
  set.seed(42)
  loo_boot <- errorest_loo_boot(x = iris_x, y = iris_y, train = MASS:::lda,
                              classify = lda_wrapper)

  # Each of the following 3 calls should result in the same error rate.
  # 1. The apparent error rate is provided, while the LOO-Boot must be computed.
  set.seed(42)
  first_error <- errorest_632plus(x = iris_x, y = iris_y, train = MASS:::lda,
                                  classify = lda_wrapper, apparent = apparent)

  # 2. The LOO-Boot error rate is provided, while the apparent must be computed.
  set.seed(42)
  second_error <- errorest_632plus(x = iris_x, y = iris_y, train = MASS:::lda,
                                   classify = lda_wrapper, loo_boot = loo_boot)

  # 3. Both error rates are provided, so the calculation is quick.
  third_error <- errorest_632plus(x = iris_x, y = iris_y, train = MASS:::lda,
                                  classify = lda_wrapper, apparent = apparent,
                                  loo_boot = loo_boot)

  # This value was computed previously.
  expected_estimate <- 0.02194472

  expect_equal(first_error, second_error)
  expect_equal(first_error, third_error)
  expect_equal(first_error, expected_estimate, tolerance = 1e-6)
})


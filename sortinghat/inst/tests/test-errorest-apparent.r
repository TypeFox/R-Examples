library('testthat')
library('sortinghat')

context("Apparent Error Rate")

test_that("Error rate works correctly on the Iris data set using MASS:::lda", {
  require('MASS')
  iris_x <- data.matrix(iris[, -5])
  iris_y <- iris[, 5]

  lda_wrapper <- function(object, newdata) {
    predict(object, newdata)$class
  }

  error_rate <- errorest_apparent(x = iris_x, y = iris_y, train = MASS:::lda,
                                  classify = lda_wrapper)

  lda_out <- MASS:::lda(x = iris_x, grouping = iris_y)
  lda_classifications <- predict(lda_out, newdata = iris_x)$class
  expected_estimate <- mean(lda_classifications != iris_y)

  expect_equal(error_rate, expected_estimate)
})


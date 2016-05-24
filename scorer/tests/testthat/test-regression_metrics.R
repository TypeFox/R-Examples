library(scorer)
context("regression metrics")

# absolute_error
# percent_error
# log_error
# squared_error
# squared_log_error
# absolute_percent_error
# mean_error
# mean_absolute_error
# median_absolute_error
# mean_percent_error
# median_percent_error
# mean_squared_error
# median_squared_error
# mean_squared_log_error
# median_squared_log_error
# mean_absolute_percent_error
# median_absolute_percent_error
# symmetric_mean_absolute_percent_error
# symmetric_median_absolute_percent_error
# mean_absolute_scaled_error
# total_variance_score
# explained_variance_score
# unexplained_variance_score
# r2_score

n <- 10000
x <- 1:n

test_that("absolute_error() produces correct output.", {
  expect_equal(absolute_error(x, x), rep(0, n))
  expect_equal(absolute_error(0, 1), 1)
  expect_equal(absolute_error(1, 0), 1)
  expect_equal(absolute_error(FALSE, TRUE), 1)
  expect_equal(absolute_error(TRUE, FALSE), 1)
})

test_that("absolute_error() produces correct output classes and types.", {
  expect_is(absolute_error(runif(n), runif(n)), "numeric")
})

test_that("absolute_error() raises correct messages, warnings, and errors.", {
  expect_error(absolute_error("a", "b"))
})

test_that("percent_error() produces correct output.", {
  expect_equal(percent_error(x, x), rep(0, n))
  expect_equal(percent_error(0, 1), -Inf)
  expect_equal(percent_error(1, 0), 1)
  expect_equal(percent_error(FALSE, TRUE), -Inf)
  expect_equal(percent_error(TRUE, FALSE), 1)
})

test_that("percent_error() produces correct output classes and types.", {
  expect_is(percent_error(runif(n), runif(n)), "numeric")
})

test_that("percent_error() raises correct messages, warnings, and errors.", {
  expect_error(percent_error("a", "b"))
})

test_that("log_error() produces correct output.", {
  expect_equal(log_error(2 * exp(1), exp(1)), 1)
  expect_equal(log_error(0, 1), NaN)
  expect_equal(log_error(1, 0), 0)
  expect_equal(log_error(FALSE, TRUE), NaN)
  expect_equal(log_error(TRUE, FALSE), 0)
})

test_that("log_error() produces correct output classes and types.", {
  expect_is(log_error(runif(n), runif(n)), "numeric")
})

test_that("log_error() raises correct messages, warnings, and errors.", {
  expect_error(log_error("a", "b"))
})

test_that("squared_error() produces correct output.", {
  expect_equal(squared_error(x, x), rep(0, n))
  expect_equal(squared_error(0, 1), 1)
  expect_equal(squared_error(1, 0), 1)
  expect_equal(squared_error(FALSE, TRUE), 1)
  expect_equal(squared_error(TRUE, FALSE), 1)
})

test_that("squared_error() produces correct output classes and types.", {
  expect_is(squared_error(runif(n), runif(n)), "numeric")
})

test_that("squared_error() raises correct messages, warnings, and errors.", {
  expect_error(squared_error("a", "b"))
})

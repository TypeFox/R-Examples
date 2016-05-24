test_that("test.is_leaf.a_false_leaf.returns_false", {
  x <- structure(list(), leaf = FALSE)
  expect_false(is_leaf(x))
})

test_that("test.is_leaf.a_leaf.returns_true", {
  x <- structure(list(), leaf = TRUE)
  expect_true(is_leaf(x))
})

test_that("test.is_leaf.a_non_logical_leaf.returns_false", {
  x <- structure(list(), leaf = 1:10)
  expect_false(is_leaf(x))
})

test_that("test.is_leaf.a_null_leaf.returns_false", {
  x <- list()
  expect_false(is_leaf(x))
})

test_that("test.is_mts.a_multivariate_time_series.returns_true", {
  expect_true(is_mts(ts(matrix(1:12, nrow = 3))))
})

test_that("test.is_mts.a_univariate_time_series.returns_false", {
  expect_false(is_mts(ts(1:12)))
})

test_that("test.is_stepfun.a_regular_function.returns_false", {
  expect_false(is_stepfun(function() {}))
})

test_that("test.is_stepfun.a_step_function.returns_true", {
  x <- stepfun(1:3, c(1, 2, 4, 3), f = 0)
  expect_true(is_stepfun(x))
})

test_that("test.is_stepfun.not_a_function.returns_false", {
  expect_false(is_stepfun(call("sin", "pi")))
})

test_that("test.is_ts.a_time_series.returns_true", {
  expect_true(is_ts(ts(1:10)))
})

test_that("test.is_ts.not_a_time_series.returns_false", {
  expect_false(is_ts(1:10))
})

test_that("test.is_tskernel.a_time_series_kernel.returns_true", {
  expect_true(is_tskernel(kernel("daniell", 10)))
})

test_that("test.is_tskernel.not_a_time_series_kernel.returns_false", {
  expect_false(is_tskernel(1:10))
})

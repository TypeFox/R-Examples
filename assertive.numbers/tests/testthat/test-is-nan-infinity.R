# TODO: tests for infiniteness

test_that("test.is_not_nan.nan.returns_false", {
  expect_false(is_not_nan(NaN))
})

test_that("test.is_not_nan.not_na.returns_true", {
  expect_true(is_not_nan(1))
})


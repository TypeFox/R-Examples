
context("Bootstrap functions")

x <- c(0.1, 1, 0.5, 0.5, 0.6, 0.71)

test_that("sampling_from_1d_gaussian_mixture", {
  xb <- rDEA:::sampling_from_1d_gaussian_mixture(3, x, 0.1)
  expect_equal( length(xb), 3L )
  expect_true( is.numeric(xb) )
})

test_that("sampling_with_reflection", {
  xb <- rDEA:::sampling_with_reflection(3, x, 0.1, 0.25, 0.5)
  expect_equal( length(xb), 3L )
  expect_true( is.numeric(xb) )
  expect_true( all(xb <= 1.0) )
})

test_that("sampling_delta_with_reflection", {
  xb <- rDEA:::sampling_delta_with_reflection(3, 1/x, 0.1, 0.25, 1.3)
  expect_equal( length(xb), 3L )
  expect_true( is.numeric(xb) )
  expect_true( all(xb >= 1.0) )
})


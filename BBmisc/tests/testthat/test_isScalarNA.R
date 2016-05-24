context("isScalarNA")

test_that("isScalarNA", {
  expect_true(isScalarNA(NA))
  expect_false(isScalarNA(1))
  expect_false(isScalarNA(iris))
  expect_false(isScalarNA(NULL))
  expect_false(isScalarNA("NA"))
})
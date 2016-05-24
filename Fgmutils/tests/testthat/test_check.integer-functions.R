context("check.integer")

test_that("check.integer return", {
  expect_equal(check.integer(5), TRUE)
  expect_equal(check.integer(5.5), FALSE)
  expect_equal(check.integer("abc"), FALSE)
})

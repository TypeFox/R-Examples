# Tests are for the various utility functions in the package
context("Utility functions")


test_that("utility functions work correctly", {

  # Expectations
  expect_identical(AnyNA(c(1, "a", NA, NULL)), TRUE)
  expect_identical(makeData(1:3, label = "a"), data.frame("a" = 1:3))

})
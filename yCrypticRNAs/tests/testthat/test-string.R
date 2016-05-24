context("Parsing and substring a character vector")

test_that("substring from the end of the string.", {
  expect_equal(substr_left("output.txt", -4), "output")
  expect_equal(substr_left("output.txt", 6), "output")
})

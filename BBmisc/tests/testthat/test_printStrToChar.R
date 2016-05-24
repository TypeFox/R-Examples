context("printStrToChar")

test_that("printStrToChar", {
  x = 1L
  s = printStrToChar(x, collapse=NULL)
  expect_equal(s, " int 1")
  s = printStrToChar(iris, collapse=NULL)
  expect_true(is.character(s) && length(s) == 6)
  s = printStrToChar(iris)
  expect_true(is.character(s) && length(s) == 1)
})
context("printToChar")

test_that("printToChar", {
  if (!interactive()) {
  z = list()
  class(z) = "foo"
  print.foo <<- function(x, ...) catf("bar")
  s = printToChar(z)
  expect_equal(s, "bar")
  print.foo <<- function(x, ...) catf("bar\nblubb")
  s = printToChar(z)
  expect_equal(s, "bar\nblubb")
  }
})
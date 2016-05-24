context("circularstring")

test_that("circularstring works", {
  # empty
  empty <- circularstring("empty")
  expect_is(empty, "character")
  expect_equal(empty, "CIRCULARSTRING EMPTY")

  # character
  str <- "CIRCULARSTRING(1 5, 6 2, 7 3)"
  a <- circularstring("CIRCULARSTRING(1 5, 6 2, 7 3)")
  expect_is(a, "character")
  expect_match(a, "CIRCULARSTRING")
  expect_equal(a, str)

  # data.frame
  str <- "CIRCULARSTRING (-116.4 45.2, -118.0 47.0)"
  df <- data.frame(lon = c(-116.4, -118), lat = c(45.2, 47))
  b <- circularstring(df, fmt = 1)
  expect_is(b, "character")
  expect_match(b, "CIRCULARSTRING")
  expect_equal(b, str)

  # matrix
  d <- "CIRCULARSTRING (-116.4 45.2, -118.0 47.0)"
  mat <- matrix(c(-116.4,-118, 45.2, 47), ncol = 2)
  d <- circularstring(mat, fmt = 1)
  expect_is(d, "character")
  expect_match(d, "CIRCULARSTRING")
  expect_equal(d, str)

  # list
  str <- "CIRCULARSTRING (1.00 5.00, 6.00 2.00, 7.00 3.00)"
  x <- list(c(1, 5), c(6, 2), c(7, 3))
  e <- circularstring(x, fmt = 2)
  expect_is(e, "character")
  expect_match(e, "CIRCULARSTRING")
  expect_equal(e, str)
})

test_that("circularstring fails correctly", {
  expect_error(circularstring(-116.4), "no applicable method")
  expect_error(circularstring(), "no applicable method")
  expect_error(circularstring(NA), "no applicable method")
  expect_error(circularstring("a", "Adf", mtcars), "The following strings are not WKT")
})

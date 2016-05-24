context("multipolygon")

test_that("multipolygon works", {
  ## empty
  empty <- multipolygon("empty")
  expect_is(empty, "character")
  expect_equal(empty, "MULTIPOLYGON EMPTY")

  ## from data.frame
  str <- "MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))"
  df <- data.frame(long = c(30, 45, 10, 30), lat = c(20, 40, 40, 20))
  df2 <- data.frame(long = c(15, 40, 10, 5, 15), lat = c(5, 10, 20, 10, 5))
  multipolydf <- multipolygon(df, df2, fmt = 0)
  expect_is(multipolydf, "character")
  expect_equal(multipolydf, str)

  ## from a matrix
  str <- "MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))"
  mat <- matrix(c(df$long, df$lat), ncol = 2)
  mat2 <- matrix(c(df2$long, df2$lat), ncol = 2)
  multipolymat <- multipolygon(mat, mat2, fmt = 0)
  expect_is(multipolymat, "character")
  expect_equal(multipolymat, str)

  ## from a list
  str <- "MULTIPOLYGON (((30 20, 45 40, 10 40, 30 20)), ((15 5, 40 10, 10 20, 5 10, 15 5)))"
  multipolylist <- multipolygon(list(c(30, 20), c(45, 40), c(10, 40), c(30, 20)),
     list(c(15, 5), c(40, 10), c(10, 20), c(5, 10), c(15, 5)), fmt = 0)
  expect_is(multipolylist, "character")
  expect_equal(multipolylist, str)
})

test_that("multipolygon fails correctly", {
  expect_error(multipolygon(-116.4), "no applicable method")
  expect_error(multipolygon(), "no applicable method")
  expect_error(multipolygon(NA), "no applicable method")
  expect_error(multipolygon("a", "Adf"), "The following strings are not WKT")
})

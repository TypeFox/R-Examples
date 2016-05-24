context("multilinestring")

test_that("multilinestring works", {
  # empty
  empty <- multipolygon("empty")
  expect_is(empty, "character")
  expect_equal(empty, "MULTIPOLYGON EMPTY")

  # character string input
  x <- "MULTILINESTRING ((30 20, 45 40, 10 40), (15 5, 40 10, 10 20))"
  multiline_char <- multilinestring(x)
  expect_is(multiline_char, "character")
  expect_match(multiline_char, "MULTILINESTRING")
  expect_equal(multiline_char, x)

  ## single point
  pt_a <- point(-116.4, 45.2)
  pt_b <- point(0, 1)
  expect_is(pt_a, "character")
  expect_is(pt_b, "character")
  expect_match(pt_a, "POINT")
  expect_match(pt_b, "POINT")
  expect_equal(point(-116.4, 45.2), "POINT (-116.4000000000000057 45.2000000000000028)")
  expect_equal(point(-116.4, 45.2, fmt = 1), "POINT (-116.4 45.2)")

  ## single point, from data.frame
  df <- data.frame(lon = -116.4, lat = 45.2)
  df_a <- point(df, fmt = 2)
  expect_is(df_a, "character")
  expect_is(df_a[1], "character")
  expect_equal(df_a[1], "POINT (-116.40 45.20)")

  ## many points, from data.frame
  df2 <- us_cities[1:5, c('lat', 'long')]
  df_b <- point(df2)
  expect_is(df_b, "character")
  expect_is(df_b[1], "character")
  expect_equal(df_b[1], "POINT (32.4500000000000028 -99.7399999999999949)")

  ## many points, from matrix
  str <- "MULTILINESTRING ((30 20, 45 40, 10 40), (15 5, 40 10, 10 20))"
  df <- data.frame(long = c(30, 45, 10), lat = c(20, 40, 40))
  df2 <- data.frame(long = c(15, 40, 10), lat = c(5, 10, 20))
  mat <- matrix(c(df$long, df$lat), ncol = 2)
  mat2 <- matrix(c(df2$long, df2$lat), ncol = 2)
  multiline_mat <- multilinestring(mat, mat2, fmt = 0)
  expect_is(multiline_mat, "character")
  expect_equal(multiline_mat, str)

  ## multilinestring, from a list
  str <- "MULTILINESTRING ((30 20, 45 40, 10 40), (15 5, 40 10, 10 20))"
  x1 <- list(c(30, 20), c(45, 40), c(10, 40))
  x2 <- list(c(15, 5), c(40, 10), c(10, 20))
  ml_list <- multilinestring(x1, x2, fmt = 0)
  expect_is(ml_list, "character")
  expect_equal(ml_list, str)

  ## test make1multilinestr helper fxn
  z <- list(c(30, 20), c(45, 40), c(10, 40))
  expect_equal(make1multilinestr(z, 0), "(30 20, 45 40, 10 40)")
})

test_that("multilinestring fails correctly", {
  expect_error(multilinestring(-116.4), "no applicable method")
  expect_error(multilinestring(), "no applicable method")
  expect_error(multilinestring(NA), "no applicable method")
  expect_error(multilinestring("a", "Adf"), "The following strings are not WKT")
})

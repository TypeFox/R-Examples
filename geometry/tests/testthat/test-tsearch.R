context("tsearch")
test_that("tsearch gives the expected output", {
  x <- c(-1, -1, 1)
  y <- c(-1, 1, -1)
  p <- cbind(x, y)
  tri <- matrix(c(1, 2, 3), 1, 3)
  ## Should be in triangle #1
  ts <- tsearch(x, y, tri, -1, -1)
  expect_that(ts, equals(1))
  ## Should be in triangle #1
  ts <- tsearch(x, y, tri, 1, -1)
  expect_that(ts, equals(1))
  ## Should be in triangle #1
  ts <- tsearch(x, y, tri, -1, 1)
  expect_that(ts, equals(1))
  ## Centroid
  ts <- tsearch(x, y, tri, -1/3, -1/3)
  expect_that(ts, equals(1))
  ## Should be outside triangle #1, so should return NA
  ts <- tsearch(x, y, tri, 1, 1)
  expect_true(is.na(ts))
})


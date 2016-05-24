# :vim set filetype=R
context("segment")
test_that("Numeric vector", {
  act <- segment(1:5)
  exp <- data.frame(a=1:4, b=2:5)
  expect_equal(act, exp)
})

test_that("Numeric vector with padding", {
  act <- segment(1:5, do.pad=TRUE)
  exp <- data.frame(a=c(NA,1:5), b=c(1:5,NA))
  expect_equal(act, exp)
})

test_that("Date vector", {
  d <- Sys.Date() + 0:4
  act <- segment(d)
  exp <- data.frame(a=d[1:4], b=d[2:5])
  expect_equal(act, exp)
})

test_that("Disallow 2D data structures", {
  x <- matrix(rnorm(10), ncol=2)
  expect_error(segment(x), "No valid function for")
})


context("item")
test_that("item works for a small sequence", {
  v <- 1:10
  expect_equal(item(v, 5), 5) 
})

test_that("item works for a small sequence with bad index.", {
  v <- 1:10
  expect_equal(item(v, 0), NA) 
})

test_that("x can not be a 2-d array", {
  x <- matrix(rnorm(10), ncol=2)
  expect_error(item(x, 1), "No valid function for")
})



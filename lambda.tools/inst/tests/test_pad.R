context("pad")
test_that("pad pads corrrect number of NAs to head", {
  x <- 1:5
  y <- c(rep(NA, 5), 1:5)
  expect_equal(pad(x, 5), y)
})

test_that("pad pads corrrect number of NAs to tail", {
  x <- 1:5
  y <- c(1:5, rep(NA, 5))
  expect_equal(pad(x, 0, 5), y)
})

test_that("pad pads corrrect number of NAs to head and tail", {
  x <- 1:5
  y <- c(rep(NA, 5), 1:5, rep(NA, 5))
  expect_equal(pad(x, 5, 5), y)
})

test_that("pad pads corrrect number of 0s to head and tail", {
  x <- 1:5
  y <- c(rep(0, 5), 1:5, rep(0, 5))
  expect_equal(pad(x, 5, 5, default=0), y)
})

test_that("x can not be a 2-d array.", {
  x <- matrix(rnorm(10), ncol=2)
  expect_error(pad(x, 10), "No valid function for")
})



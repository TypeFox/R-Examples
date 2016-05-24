# :vim set filetype=R
context("confine")
test_that("Basic behavior is correct", {
  x <- -2.5 : 2.5
  y <- c(-1, -1, -.5, .5, 1, 1)
  expect_equal(confine(x), y)
})

test_that("The min.level and max.level options work", {
  x <- -4.5 : 4.5
  y <- c(-2, -2, -2, -1.5, -.5, .5, 1.5, 2, 2, 2)
  expect_equal(confine(x, max.level=2, min.level=-2), y)
})

test_that("Confine fails when min.level >= max.level", {
  x <- 1:100
  y <- rnorm(100, sd=4)

  expect_error(confine(y, min.len=1, max.len=-1), "No valid function for")
  expect_error(confine(y, min.len=1, max.len=1), "No valid function for")
})

test_that("Confine fails for character input", {
  x <- 'abcdefghijk'

  expect_error(confine(x), "No valid function for")
  expect_error(confine(x, min.len='b', max.len='h'), "No valid function for")
})



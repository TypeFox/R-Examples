context("is.scalar")
test_that("scalar argument returns true", {
  expect_true(is.scalar(10))
})

test_that("Non-scalar argument returns false.", {
  expect_true(!is.scalar(c(10, 10)))
})


context("onlyif")
test_that("onlyif works with positive condition", {
  x <- 1:5
  y <- onlyif(length(x) < 10, function(y) pad(y, 10 - length(y)), x)
  expect_equal(y, c(NA, NA, NA, NA, NA, 1, 2, 3, 4, 5))
})

test_that("onlyif works with negative condition", {
  x <- 1:20
  y <- onlyif(length(x) < 10, function(y) fold(x, function(x, y) x+y), x)
  expect_equal(y, x)
})

test_that("fn argument must be a function.", {
  x <- rnorm(5)
  constant <- 'string' 
  expect_equal(onlyif(TRUE, constant, x), "string")
})


context("use_default")
test_that("use_default works", {
  x <- c(1, 2, 3, NA, NA)
  y <- map(x, function(a) use_default(a, 0))
  expect_equal(y, c(1, 2, 3, 0, 0))
})

test_that("Non-scalar arguement throws an error.", {
  
  expect_error(use_default(c(1, 1), 0), "No valid function for")
})

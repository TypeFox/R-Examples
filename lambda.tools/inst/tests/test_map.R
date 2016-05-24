# :vim set filetype=R
context("1D map")
test_that("Basic functionality", {
  x <- 1:10
  act <- map(x, function(a) a^2)
  exp <- x^2
  expect_equal(act, exp)
})

test_that("Custom accumulator", {
  x <- 1:10
  act <- map(x, function(a) a^2, acc=list())
  exp <- as.list(x^2)
  expect_equal(act, exp)
})

test_that("Invalid functions are not allowed", {
  x <- 1:10
  f <- 'I am not a function.'
  expect_error(map(x, f))
})

test_that("NAs okay in input sequence", {
  x <- c(-5:5, NA, NA, 7:10) 
  act <- map(x, function(a) abs(a))
  exp <- abs(x)
  expect_equal(act, exp)
})


context("2D map")
test_that("Basic functionality", {
  x <- matrix(1:10, ncol=2)
  act <- map(x, function(a) sum(a))
  exp <- apply(x, 2, sum)
  expect_equal(act, exp)
})

test_that("Use a data.frame as input", {
  x <- data.frame(a=1:10, b=11:20)
  act <- map(x, function(a) sum(a))
  exp <- c(55, 155)
  expect_equal(act, exp)
})

test_that("Use a data.frame as input and accumulate in a list", {
  x <- data.frame(a=1:10, b=11:20)
  act <- map(x, function(a) sum(a), acc=list())
  exp <- list(55, 155)
  expect_equal(act, exp)
})


context("1D maprange")
test_that("maprange with window size a multiple of length(x)",{
  x <- 1:10
  y <- maprange(x, 2, function(a) sum(a))
  expect_equal(y, c(3, 5, 7, 9, 11, 13, 15, 17, 19))
})

test_that("maprange with window size not a multiple of length(x)",{
  x <- 1:10
  y <- maprange(x, 3, function(a) sum(a))
  expect_equal(y, c(6, 9, 12, 15, 18, 21, 24, 27))
})

test_that("maprange with window size not a multiple of length(x) and do.pad",{
  x <- 1:10
  y <- maprange(x, 3, function(a) sum(a), TRUE)
  expect_equal(y, c(NA, NA, 6, 9, 12, 15, 18, 21, 24, 27))
})

test_that("maprange with x a vector and with window size > length(x)",{
  x <- 1:10
  expect_error(maprange(x, 11, function(a) sum(a), "No valid function for"))
})


context("2D maprange")
# TODO: Add tests here


context("1D mapblock")
test_that("Basic functionality",{
  x <- 1:10
  act <- mapblock(x, 2, mean)
  exp <- apply(matrix(x,nrow=2), 2, mean)
  expect_equal(act, exp)
})

test_that("Window does not divide length of x",{
  x <- 1:10
  act <- mapblock(x, 3, mean)
  exp <- c(apply(matrix(x[1:9],nrow=3), 2, mean), 10)
  expect_equal(act, exp)
})

test_that("Window longer than length of x",{
  x <- 1:10
  act <- mapblock(x, 11, mean)
  exp <- mean(x)
  expect_equal(act, exp)
})


context("2D mapblock")
test_that("mapblock with x a matrix block size of one",{
  m <- matrix(1:12, ncol=2)
  y <- mapblock(m, 1, sum)
  expect_equal(y, c(21, 57))
})

test_that("mapblock with x a matrix block size equal to ncol(m)",{
  m <- matrix(1:12, ncol=2)
  y <- mapblock(m, 2, sum)
  expect_equal(y, c(78)) 
})

test_that("mapblock with x a data.frame block size of one",{
  df <- data.frame(col1=1:6, col2=7:12)
  y <- mapblock(df, 1, sum)
  expect_equal(y, c(21, 57))
})



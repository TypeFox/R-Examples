# :vim set filetype=R
context("slice")
test_that("slice a sequence with inclusive flag set to FALSE", {
  act <- slice(1:10, 5)
  exp <- list(c(1:5), c(6:10)) 
  expect_equal(act, exp)
})

test_that("slice a sequence with inclusive flag set to TRUE", {
  act <- slice(1:10, 5, TRUE)
  exp <- list(c(1:5), c(5:10)) 
  expect_equal(act, exp)
})

test_that("slice a sequence and an expression", {
  x <- 1:10
  y <- list(c(1:4), c(5:10))
  expect_equal(slice(x, x < 5), y)
})

test_that("slice a sequence and another expression", {
  x <- 1:10
  y <- list(c(6:10), c(1:5))
  expect_equal(slice(x, x > 5), y)
})

test_that("slice a matrix", {
  A <- matrix(1:10, ncol=2)
  y <- list(matrix(c(1,2,3,6,7,8),ncol=2), matrix(c(3,4,5,8,9,10), ncol=2))
  expect_equal(slice(A, 3, TRUE), y)
})

test_that("slice a data.frame inclusive", {
  df <- data.frame(col1=1:10, col2=1:10)
  y <- list(df[1:5,], df[5:10,])
  expect_equal(slice(df, 5, TRUE), y)
})

test_that("slice a data.frame non-inclusive", {
  df <- data.frame(col1=1:10, col2=1:10)
  y <- list(df[1:5,], df[6:10,])
  expect_equal(slice(df, 5, FALSE), y)
})

test_that("slice a data.frame with an expression", {
  df <- data.frame(col1=1:10, col2=11:20)
  act <- slice(df, df$col1 < 5)
  exp <- list(df[1:4,], df[5:10,])
  expect_equal(act, exp)
})

test_that("pivot is larger than length of x", {
  x <- 1:50
  df <- data.frame(x=x, y=x)
  expect_error(slice(1:50, 51), "No valid function for")
  expect_error(slice(df, 51), "No valid function for")
})

test_that("expression length is not equal to length of x", {
  x <- 1:50
  x2 <- 1:100
  df <- data.frame(col1=x, col2=x) 
  df2 <- data.frame(col1=x2, col2=x2)
  expect_error(slice(x, x2 < 25), "No valid function for")
  expect_error(slice(x, x2 < 25 & x2 > 50), "No valid function for")
  expect_error(slice(df, df2$col1 < 25 & df2$col1 > 50), "No valid function for")
  expect_error(slice(df, df2[,1] < 25 & df2[,1] > 50), "No valid function for")
})



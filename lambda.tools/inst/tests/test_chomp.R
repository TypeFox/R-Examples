# :vim set filetype=R
context("chomp")
test_that("chomp works with a small sequence and default head/tail param", {
  x <- 1:10
  y <- 2:9
  expect_equal(chomp(x), y)
})

test_that("chomp works with a small sequence with different head/tail params", {
  x <- 1:10
  y <- 3:8
  expect_equal(chomp(x, head=2, tail=2), y)
})

test_that("chomp works with a matrix", {
  m <- matrix(1:10, ncol=2)
  y <- matrix(c(2, 3, 4, 7, 8, 9), ncol=2) 
  expect_equal(chomp(m), y)
})

test_that("chomp works on a data.frame", {
  df <- data.frame(x=1:10, y=1:10)
  y <- df[3:8,]
  expect_equal(chomp(df, head=2, tail=2), y)
})

test_that("Chomp with head and tail < 0 fails.", {
  x <- 1:50
  df <- data.frame(col1=x, col2=x)
  expect_error(chomp(x, head=-1, tail=-1), "No valid function for")
  expect_error(chomp(df, head=-1, tail=-1), "No valid function for")
})

test_that("Head and tail parameters can not overlap.", {
  x <- 1:50
  df <- data.frame(col1=x, col2=x)
  expect_error(chomp(x, head=26, tail=26), "No valid function for")
  expect_error(chomp(df, head=1, tail=50), "No valid function for")
})

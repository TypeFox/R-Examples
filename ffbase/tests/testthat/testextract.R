library(testthat)

context("extract")

test_that("Extracting with logical vector works", {
  x <- 1:2
  y <- c(TRUE, FALSE)
  x_ff <- ff(x)
  y_ff <- ff(y)
  
  expect_equal(x[y], x_ff[y_ff][])
})
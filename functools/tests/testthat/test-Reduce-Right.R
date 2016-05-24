library(functools)
context("Reduce_Right()")

divide_x_by_y <- function(x, y) return(x / y)

z <- c(128, 2, 2, 2, 2, 2, 2, 2)

test_that("Produces the correct output.", {
  expect_equal(Reduce_Right(divide_x_by_y, z),
               Reduce(divide_x_by_y, z, right = TRUE))
})

test_that("Produces the correct output type.", {
  expect_is(Reduce_Right(divide_x_by_y, z), "numeric")
})

test_that("Produces the correct errors.", {
  expect_error(Reduce_Right(divide_x_by_y, z)())
})

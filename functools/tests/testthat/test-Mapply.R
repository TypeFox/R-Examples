library(functools)
context("Mapply()")

test_that("Produces the correct output.", {
  expect_equal(Mapply(1:4, 4:1, .f = rep),
               mapply(rep, 1:4, 4:1))
  expect_equal(Mapply(times = 1:4, x = 4:1, .f = rep),
               mapply(rep, times = 1:4, x = 4:1))
  expect_equal(Mapply(times = 1:4, more_args = list(x = 42), .f = rep),
               mapply(rep, times = 1:4, MoreArgs = list(x = 42)))
})

test_that("Produces the correct output type.", {
  expect_is(Mapply(1:4, 4:1, .f = rep), "list")
})

test_that("Produces the correct errors.", {
  expect_error(Mapply(rep, 1:4, 1:4))
})

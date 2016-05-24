library(functools)
context("Best()")

test_that("Produces the correct output.", {
  expect_equal(Best(1:10, function(x, y) return(x > y)), 10)
  expect_equal(Best(1:10, function(x, y) return(x < y)), 1)
  expect_equal(Best(letters, function(x, y) return(x[1] == "l")), "l")
})

test_that("Produces the correct output type.", {
  expect_is(Best(1:10, function(x, y) return(x > y)), "integer")
  expect_is(Best(1:10, function(x, y) return(x < y)), "integer")
  expect_is(Best(letters, function(x, y) return(x[1] == "l")), "character")
})

test_that("Produces the correct errors.", {
  expect_error(Best(1:10, function(x, y, z) return(y < z)))
})


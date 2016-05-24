library(functools)
context("Compact()")

foo <- list(NA_integer_, TRUE, NA_character_, 1, NULL,
            "a", NA_complex_, 4.5, FALSE, NA_real_)
bar <- list(TRUE, 1, "a", 4.5, FALSE)
test_that("Produces the correct output.", {
  expect_equal(Compact(foo), bar)
  expect_equal(Compact(bar), bar)
  expect_equal(Compact(mtcars), mtcars)
  expect_equal(Compact(1L), 1L)
  expect_equal(Compact(2.6), 2.6)
  expect_equal(Compact(TRUE), TRUE)
  expect_equal(Compact("a"), "a")
  expect_equal(Compact(NULL), NULL)
  expect_equal(Compact(NA), logical(0))
})

test_that("Produces the correct output type.", {
  expect_is(Compact(bar), "list")
  expect_is(Compact(foo), "list")
  expect_is(Compact(mtcars), "data.frame")
  expect_is(Compact(1L), "integer")
  expect_is(Compact(2.6), "numeric")
  expect_is(Compact(TRUE), "logical")
  expect_is(Compact("a"), "character")
  expect_is(Compact(NULL), "NULL")
  expect_is(Compact(NA), "logical")
})

test_that("Produces the correct errors.", {
  expect_equal(1, 1)
})


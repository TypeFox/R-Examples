library(functools)
context("Truthy()")

foo <- list(NA_integer_, TRUE, NA_character_, 1, NULL,
            "a", NA_complex_, 4.5, FALSE, NA_real_)
bar <- c(FALSE, TRUE, FALSE, FALSE, FALSE,
         FALSE, FALSE, FALSE, FALSE, FALSE)

test_that("Produces the correct output.", {
  expect_equal(unlist(lapply(foo, Truthy)), bar)
  expect_equal(Truthy(NULL), FALSE)
  expect_equal(Truthy(NA), FALSE)
  expect_warning(Truthy(Truthy), "is.na() applied to non-(list or vector) of type 'closure'", fixed = TRUE)
  expect_equal(Truthy(1), FALSE)
  expect_equal(Truthy(0), FALSE)
  expect_equal(Truthy(2.4), FALSE)
  expect_equal(Truthy("hello"), FALSE)
  expect_equal(Truthy(FALSE), FALSE)
  expect_equal(Truthy(TRUE), TRUE)
})

test_that("Produces the correct output type.", {
  expect_is(Truthy(NULL), "logical")
  expect_is(Truthy(NA), "logical")
  expect_is(Truthy(1), "logical")
  expect_is(Truthy(0), "logical")
  expect_is(Truthy(2.4), "logical")
  expect_is(Truthy("hello"), "logical")
  expect_is(Truthy(FALSE), "logical")
  expect_is(Truthy(TRUE), "logical")
})

test_that("Produces the correct errors.", {
  expect_warning(Truthy(Truthy), "is.na() applied to non-(list or vector) of type 'closure'", fixed = TRUE)
})

library(functools)
context("Always()")

test_that("Produces the correct output.", {
  expect_equal(Always(TRUE)(), TRUE)
  expect_equal(Always(FALSE)(), FALSE)
  expect_equal(Always("a")(), "a")
  expect_equal(Always(1L)(), 1L)
  expect_equal(Always(2.6)(), 2.6)
  expect_equal(Always(NULL)(), NULL)
  expect_equal(Always(NA)(), NA)
})

test_that("Produces the correct output type.", {
  expect_is(Always(), "function")
  expect_is(Always(TRUE), "function")
  expect_is(Always(TRUE)(), "logical")
})

test_that("Produces the correct errors.", {
  expect_error(Always()())
})


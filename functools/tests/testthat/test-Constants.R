library(functools)
context("Constants functions")

test_that("Produces the correct output.", {
  expect_equal(True(), TRUE)
  expect_equal(False(), FALSE)
  expect_equal(Identity(Identity), Identity)
  expect_equal(Na(), NA)
  expect_equal(Null(), NULL)
})

test_that("Produces the correct output type.", {
  expect_is(True(), "logical")
  expect_is(False(), "logical")
  expect_is(Identity(Identity), "function")
  expect_is(Na(), "logical")
  expect_is(Null(), "NULL")
})

test_that("Produces the correct errors.", {
  expect_error(True(1))
  expect_error(False(1))
  expect_error(Na(1))
  expect_error(Null(1))
})

library(functools)
context("Any()")

test_that("Produces the correct output.", {
  expect_equal(Any(mtcars, is.numeric), TRUE)
  expect_equal(Any(mtcars, is.character), FALSE)
  expect_equal(Any(list(NA, "3", NULL), is.numeric), FALSE)
  expect_equal(Any(list(NA, "3", NULL, 5), is.numeric), TRUE)
  expect_equal(Any(list(NA, 3, NULL), is.numeric), TRUE)
  expect_equal(Any(list(NA, FALSE), Identity), NA)
  expect_equal(Any(list(NA, FALSE), Identity, na.rm = TRUE), FALSE)
})

test_that("Produces the correct output type.", {
  expect_is(Any(mtcars, is.numeric), "logical")
  expect_is(Any(mtcars, is.character), "logical")
  expect_is(Any(list(NA, "3", NULL), is.numeric), "logical")
  expect_is(Any(list(NA, "3", NULL, 5), is.numeric), "logical")
  expect_is(Any(list(NA, 3, NULL), is.numeric), "logical")
  expect_is(Any(list(NA, FALSE), Identity), "logical")
  expect_is(Any(list(NA, FALSE), Identity, na.rm = TRUE), "logical")
})

test_that("Produces the correct errors.", {
  expect_error(Any(mtcars, mean), "logical")
})

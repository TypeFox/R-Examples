library(functools)
context("All()")

test_that("Produces the correct output.", {
  expect_equal(All(mtcars, is.numeric), TRUE)
  expect_equal(All(mtcars, is.character), FALSE)
  expect_equal(All(list(NA, "3", NULL), is.numeric), FALSE)
  expect_equal(All(list(NA, "3", NULL, 5), is.numeric), FALSE)
  expect_equal(All(list(NA, 3, NULL), is.numeric), FALSE)
  expect_equal(All(list(NA, TRUE), Identity), NA)
  expect_equal(All(list(NA, TRUE), Identity, na.rm = TRUE), TRUE)
})

test_that("Produces the correct output type.", {
  expect_is(All(mtcars, is.numeric), "logical")
  expect_is(All(mtcars, is.character), "logical")
  expect_is(All(list(NA, "3", NULL), is.numeric), "logical")
  expect_is(All(list(NA, "3", NULL, 5), is.numeric), "logical")
  expect_is(All(list(NA, 3, NULL), is.numeric), "logical")
  expect_is(All(list(NA, TRUE), Identity), "logical")
  expect_is(All(list(NA, TRUE), Identity, na.rm = TRUE), "logical")
})

test_that("Produces the correct errors.", {
  expect_error(All(mtcars, mean), "logical")
})


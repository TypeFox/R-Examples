context("Basic tests")

test_that("logical tests act as expected", {
  expect_that(TRUE, is_true())
  expect_that(FALSE, is_false())
})

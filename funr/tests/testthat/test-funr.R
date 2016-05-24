

library(testthat)

test_that("check return of fetch_pipes", {
  shows_message(funr())
  expect_error(funr("rnorm"), '')
})

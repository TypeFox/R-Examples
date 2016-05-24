context("topics")

test_that("fails well with no input", {
	skip_on_cran()

  expect_error(topic(), "argument \"id\" is missing")
})

test_that("fails well with non-existent page", {
	skip_on_cran()

  expect_error(topic("asfafsfadfasdfd"),
               "404 - Oops! That page doesn’t exist or is private.")
})

test_that("httr curl options work", {
	skip_on_cran()

  library("httr")
  expect_error(topic("asdfadf", config = timeout(seconds = 0.001)))
})

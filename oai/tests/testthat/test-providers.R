context("providers")

test_that("update_providers", {
  skip_on_cran()

  library("httr")

  # call would take too long, mocking
  expect_error(update_providers(config = timeout(0.001)),
               "Timeout was reached")
})

test_that("load_providers", {
  # delete providers
  stuff <- providers
  expect_is(stuff, "data.frame")
  rm(stuff)
  expect_error(get("stuff"))
})

test_that("providers data set", {
  expect_is(providers, "data.frame")
  expect_less_than(2500, NROW(providers))
})

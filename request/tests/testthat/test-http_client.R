context("http_client")

test_that("http_client works", {
  skip_on_cran()

  aa <- api("http://api.plos.org/search") %>%
    api_query(q = ecology, wt = json, fl = 'id,journal') %>%
    http_client()

  expect_is(aa, "RequestIterator")
  expect_is(aa$body(), "response")
  expect_is(aa$count(), "integer")
  expect_is(aa$parse(), "list")
})

test_that("http_client fails well", {
  skip_on_cran()

  expect_error(http_client(), "argument \"req\" is missing")
})

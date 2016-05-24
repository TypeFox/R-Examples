context("os_statemetasearch")

test_that("returns the correct", {
  skip_on_cran()

  a <- os_statemetasearch()
  b <- os_statemetasearch(state = 'ca')

  # classes
  expect_is(a, "data.frame")
  expect_is(b, "list")
  expect_is(b$ca, "list")
  expect_is(b$ca$terms, "data.frame")

  # values
  expect_named(b, "ca")
})

test_that("vectorizing works", {
  skip_on_cran()

  res <- os_statemetasearch(c('tx','nv'))
  expect_is(res, "list")
  expect_equal(length(res), 2)
  expect_named(res, c('tx', 'nv'))
})

test_that("curl options work", {
  skip_on_cran()

  expect_error(os_statemetasearch(config = httr::timeout(0.001)), "Timeout")
})

test_that("fails well", {
  skip_on_cran()

  expect_error(os_statemetasearch(4), "Not Found")
  expect_error(os_statemetasearch("adfad"), "Not Found")
})

context("test idig_meta_fields")

test_that("records field list returns", {
  testthat::skip_on_cran()
  f <- idig_meta_fields()

  expect_that(f, is_a("list"))
  expect_that(f[["data"]], is_a("list"))
  expect_that(f[["uuid"]][["type"]] == "string", is_true())
})

test_that("records indexed subset returns", {
  testthat::skip_on_cran()
  f <- idig_meta_fields(subset="indexed")
  
  expect_that(f[["data"]], is_null())
  expect_that(f[["uuid"]][["type"]] == "string", is_true())
})

test_that("records raw subset returns", {
  testthat::skip_on_cran()
  f <- idig_meta_fields(subset="raw")
  
  expect_that(f[["uuid"]], is_null())
  expect_that(f[["dwc:occurrenceID"]][["type"]] == "string", is_true())
})

test_that("media list returns", {
  testthat::skip_on_cran()
  f <- idig_meta_fields(type="media")

  expect_that(f, is_a("list"))
  expect_that(f[["data"]], is_a("list"))
  expect_that(f[["uuid"]][["type"]] == "string", is_true())
})

test_that("media indexed subset returns", {
  testthat::skip_on_cran()
  f <- idig_meta_fields(type="media", subset="indexed")
  
  expect_that(f[["data"]], is_null())
  expect_that(f[["uuid"]][["type"]] == "string", is_true())
})

test_that("media raw subset returns", {
  testthat::skip_on_cran()
  f <- idig_meta_fields(type="media", subset="raw")
  
  expect_that(f[["uuid"]], is_null())
  expect_that(f[["ac:accessURI"]][["type"]] == "string", is_true())
})
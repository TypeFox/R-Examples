context("test idig_count_records")

recordset <- "7450a9e3-ef95-4f9e-8260-09b498d2c5e6"

test_that("all records is a large number", {
  testthat::skip_on_cran()
  num <- idig_count_records()

  expect_that(num, is_a("integer"))
  expect_that((num - num) == 0, is_true())
  expect_that(num > 20 * 1000 * 1000, is_true())
})

test_that("rq searches on the endpoint is a small number", {
  testthat::skip_on_cran()
  num <- idig_count_records(rq=list("recordset"=recordset))
  
  expect_that(num, is_a("integer"))
  expect_that(num > 400 * 1000, is_true())
})
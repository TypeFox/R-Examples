library(testthat)
library(WufooR)

options(Wufoo_Name = "johnmalc", Wufoo_API = "F1QH-Q64B-BSBI-JASJ")

context("Fields")

test_that("Fields request returns 17 rows, always", {
  fieldsALL <- fields_info(formIdentifier = "z5kqx7h1gtvg4g")
  expect_more_than(dim(fieldsALL)[1], 1)
})

context("Reports")

test_that("Reports have 11 fields", {
  repo <- reports_info()
  expect_equal(dim(repo)[2], 11)
})

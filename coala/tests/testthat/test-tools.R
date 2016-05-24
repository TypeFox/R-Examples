context("Tools")

test_that("it checks for packages", {
  skip_on_cran()
  expect_true(require_package("coala"))
  expect_error(require_package("2l3ihjrpaiwhf"))
})

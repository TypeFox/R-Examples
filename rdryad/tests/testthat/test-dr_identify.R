context("oai-pmh - dr_identify")

test_that("dr_identify works", {
  skip_on_cran()

  aa <- dr_identify()

  expect_is(aa, "data.frame")
  expect_is(aa$repositoryName, "character")
  expect_true(grepl("Dryad", aa$repositoryName))
})


test_that("dr_identify curl options work", {
  skip_on_cran()

  library("httr")
  expect_error(dr_identify(config = timeout(0.001)), "Timeout was reached")
})

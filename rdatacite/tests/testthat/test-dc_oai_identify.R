context("dc_oai_identify")

test_that("dc_oai_identify - default uses datacite", {
  skip_on_cran()

  aa <- dc_oai_identify()

  expect_is(aa, "data.frame")
  expect_match(aa$repositoryName, "DataCite")
  expect_match(aa$baseURL, "oai.datacite.org")
})

test_that("dc_oai_identify - httr options work", {
  skip_on_cran()

  library("httr")
  expect_error(dc_oai_identify(config = timeout(0.01)))
})

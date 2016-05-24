context("dc_oai_listsets")

test_that("dc_oai_listsets", {
  skip_on_cran()

  aa <- dc_oai_listsets()

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_is(aa$setSpec, "character")
  expect_is(aa$setName, "character")
})

test_that("dc_oai_listsets - curl options", {
  skip_on_cran()

  library("httr")

  expect_error(dc_oai_listsets(config = timeout(0.001)), "Timeout was reached")
})

test_that("dc_oai_listsets fails well", {
  skip_on_cran()

  expect_error(dc_oai_listsets(token = 454),
               "The value of the resumptionToken argument is invalid or expired")
  expect_error(dc_oai_listsets("stuff"),
               "The value of the resumptionToken argument is invalid or expired")
})

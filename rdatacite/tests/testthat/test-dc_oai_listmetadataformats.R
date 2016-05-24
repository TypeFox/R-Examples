context("dc_oai_listmetadataformats")

test_that("dc_oai_listmetadataformats - basic functionality works", {
  skip_on_cran()

  aa <- dc_oai_listmetadataformats()

  expect_is(aa, "data.frame")
  expect_is(aa$metadataPrefix, "character")
  expect_is(aa$schema, "character")
  expect_is(aa$metadataNamespace, "character")

  expect_named(aa, c('metadataPrefix', 'schema', 'metadataNamespace'))
})

test_that("dc_oai_listmetadataformats - no formats avail. vs. avail", {
  skip_on_cran()

  aa <- dc_oai_listmetadataformats(id = "oai:oai.datacite.org:22")
  bb <- dc_oai_listmetadataformats(id = "oai:oai.datacite.org:6718729")

  expect_null(aa[[1]])
  expect_is(bb[[1]], "data.frame")
})

test_that("dc_oai_listmetadataformats - curl options", {
  skip_on_cran()

  library("httr")

  expect_error(dc_oai_listmetadataformats(config = timeout(0.001)), "Timeout was reached")
})

test_that("dc_oai_listmetadataformats fails well", {
  skip_on_cran()

  expect_null(dc_oai_listmetadataformats(id = "adfadfsdf")[[1]])
})

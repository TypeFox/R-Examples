context("oai-pmh - dr_list_metadata_formats")

test_that("dr_list_metadata_formats works", {
  skip_on_cran()

  aa <- dr_list_metadata_formats()

  expect_is(aa, "data.frame")
  expect_is(aa$metadataPrefix, "character")
  expect_true(any(grepl("rdf", aa$metadataPrefix)))
})

test_that("dr_list_metadata_formats curl options work", {
  skip_on_cran()

  library("httr")
  expect_error(dr_identify(config = timeout(0.001)), "Timeout was reached")
})

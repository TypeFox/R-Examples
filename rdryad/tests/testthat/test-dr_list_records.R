context("oai-pmh - dr_list_records")

test_that("dr_list_records works - output formats", {
  skip_on_cran()

  aa <- dr_list_records(from = '2010-01-01', until = '2010-09-10')
  bb <- dr_list_records(from = '2010-01-01', until = '2010-09-10', as = "list")
  cc <- dr_list_records(from = '2010-01-01', until = '2010-09-10', as = "raw")

  expect_is(aa, "data.frame")
  expect_is(bb, "list")
  expect_is(cc, "character")

  expect_named(aa, c('identifier', 'datestamp', 'setSpec'))

  library("xml2")
  expect_is(xml2::read_xml(cc), "xml_document")

  expect_named(bb[[1]], c('headers', 'metadata'))
})

test_that("dr_list_records fails well", {
  skip_on_cran()

  expect_error(dr_list_records(prefix = 5), "\"5\" is not supported")
  expect_error(dr_list_records(from = "the"), "The request includes illegal arguments")
  expect_error(dr_list_records(until = "adfafdfd"), "The request includes illegal arguments")
  expect_error(dr_list_records(set = 344), "The request includes illegal arguments")
  # expect_error(dr_list_records(as = 5), "The request includes illegal arguments") FIXME
})

test_that("dr_list_records curl options work", {
  skip_on_cran()

  library("httr")
  expect_error(dr_list_records(config = timeout(0.001)), "Timeout was reached")
})

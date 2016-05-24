context("oai-pmh - dr_get_records")

id <- 'oai:datadryad.org:10255/dryad.8820'
handles <- c('10255/dryad.36217', '10255/dryad.86943', '10255/dryad.84720', '10255/dryad.34100')
ids <- paste0('oai:datadryad.org:', handles)

test_that("dr_get_records works", {
  skip_on_cran()

  aa <- dr_get_records(id)
  bb <- dr_get_records(ids)

  expect_is(aa, "data.frame")
  expect_is(aa$identifier, "character")
  expect_equal(NROW(aa), length(id))
  expect_true(grepl("dryad", aa$identifier))

  expect_is(bb, "data.frame")
  expect_is(bb$identifier, "character")
  expect_equal(NROW(bb), length(ids))
  expect_true(any(grepl("dryad", bb$identifier)))
})

test_that("dr_get_records fails as expected", {
  skip_on_cran()

  expect_error(dr_get_records(), "argument \"ids\" is missing, with no default")
  expect_error(dr_get_records(id, prefix = 44),
               "\"44\" is not supported")

  library("httr")
  expect_error(dr_get_records(id, config = timeout(0.001)), "Timeout was reached")
})

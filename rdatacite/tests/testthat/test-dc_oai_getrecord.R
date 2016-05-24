context("dc_oai_getrecord")

test_that("dc_oai_getrecord - basic functionality works", {
  skip_on_cran()

  aa <- dc_oai_getrecord("oai:oai.datacite.org:32255")

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_is(aa$identifier, "character")
  expect_is(aa$title, "character")
})

test_that("dc_oai_getrecord - many record Ids input works", {
  skip_on_cran()

  recs <- c("oai:oai.datacite.org:32255", "oai:oai.datacite.org:32325")
  aa <- dc_oai_getrecord(recs)

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_is(aa$identifier, "character")
  expect_is(aa$title, "character")
  expect_equal(NROW(aa), 2)
  expect_equal(aa$identifier, recs)
})

test_that("dc_oai_getrecord fails well", {
  skip_on_cran()

  expect_error(dc_oai_getrecord(),
               "argument \"id\" is missing, with no default")
  expect_error(dc_oai_getrecord('5000000000000asfaffs'),
               "is unknown or illegal in this repository")
})

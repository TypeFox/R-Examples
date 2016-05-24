context("dc_oai_listrecords")

test_that("dc_oai_listrecords works", {
  skip_on_cran()

  a <- dc_oai_listrecords(from = '2011-06-01T', until = '2011-07-01T')

  expect_is(a, "data.frame")
  expect_is(a, "oai_df")
  expect_is(a$identifier, "character")
})

test_that("dc_oai_listrecords fails well", {
  skip_on_cran()

  expect_error(dc_oai_listrecords(from = '2011-06-01T', until = 'adffdsadsf'),
               "The request includes illegal arguments")
  expect_error(dc_oai_listrecords(from = '2011-06-01T', until = 5),
               "The request includes illegal arguments")
  expect_error(dc_oai_listrecords(url = 5), "One or more of your URLs")
})

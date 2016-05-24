context("id")

test_that("id - default uses datacite", {
  skip_on_cran()

  aa <- id("http://oai.datacite.org/oai")

  expect_is(aa, "data.frame")
  expect_match(aa$repositoryName, "DataCite")
  expect_match(aa$baseURL, "oai.datacite.org")
})

test_that("id - url param works", {
  skip_on_cran()

  aa <- id("http://export.arxiv.org/oai2")
  bb <- id("http://pub.bsalut.net/cgi/oai2.cgi")
  cc <- id("http://www.diva-portal.org/oai/OAI")

  expect_is(aa, "data.frame")
  expect_is(bb, "data.frame")
  expect_is(cc, "data.frame")

  expect_equal(NROW(aa), 1)
  expect_equal(NROW(bb), 1)
  expect_equal(NROW(cc), 1)
})

test_that("id fails well", {
  skip_on_cran()

  expect_error(id(),
               "argument \"url\" is missing")
  expect_error(id("things"),
               "One or more of your URLs")
})

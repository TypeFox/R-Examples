context("count_identifiers")

test_that("count_identifiers - basic functionality", {
  skip_on_cran()

  aa <- count_identifiers()

  expect_is(aa, "data.frame")
  expect_is(aa$url, "character")
  expect_is(aa$count, "numeric")
  expect_match(aa$url, "oai.datacite.org/oai")
})

test_that("count_identifiers - works with many input urls", {
  skip_on_cran()

  aa <- count_identifiers(c(
    "http://oai.datacite.org/oai",
    "http://archivesic.ccsd.cnrs.fr/oai/oai.php"
  ))

  expect_is(aa, "data.frame")
  expect_is(aa$url, "character")
  expect_is(aa$count, "numeric")

  expect_equal(NROW(aa), 2)
})

test_that("count_identifiers - prefix param works", {
  skip_on_cran()

  aa <- count_identifiers("http://oai.datacite.org/oai", prefix = "oai_dc")
  bb <- count_identifiers("http://oai.datacite.org/oai", prefix = "oai_datacite")

  expect_is(aa, "data.frame")
  expect_is(bb, "data.frame")

  expect_equal(aa, bb)
})

test_that("count_identifiers - curl options", {
  skip_on_cran()

  library("httr")

  expect_error(count_identifiers(config = timeout(0.001)), "Timeout was reached")
})

test_that("count_identifiers fails well", {
  skip_on_cran()

  expect_error(list_sets(token = 454),
               "The value of the resumptionToken argument is invalid or expired")
  expect_error(list_sets("stuff"),
               "One or more of your URLs")
})

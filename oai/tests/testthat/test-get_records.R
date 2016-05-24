context("get_records")

test_that("get_records - basic functionality works", {
  skip_on_cran()

  aa <- get_records("oai:oai.datacite.org:32255")

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_is(aa$identifier, "character")
  expect_is(aa$title, "character")
})

test_that("get_records - many record Ids input works", {
  skip_on_cran()

  recs <- c("oai:oai.datacite.org:32255", "oai:oai.datacite.org:32325")
  aa <- get_records(recs)

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_is(aa$identifier, "character")
  expect_is(aa$title, "character")
  expect_equal(NROW(aa), 2)
  expect_equal(aa$identifier, recs)
})

test_that("get_records fails well", {
  skip_on_cran()

  expect_error(get_records(),
               "argument \"ids\" is missing, with no default")
  expect_error(get_records(5),
               "\"5\" is unknown or illegal in this repository")
  expect_error(get_records("oai:oai.datacite.org:32255", url = "stuff"),
               "One or more of your URLs")
})

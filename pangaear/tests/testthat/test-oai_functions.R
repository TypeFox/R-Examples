context("oai functions")

test_that("pg_identify() works", {
  skip_on_cran()

  aa <- pg_identify()

  expect_is(aa, "pg_identify")
  expect_is(aa$repositoryName, "character")
})

test_that("pg_list_sets() works", {
  skip_on_cran()

  aa <- pg_list_sets()

  expect_is(aa, "oai_df")
  expect_is(aa$setSpec, "character")
  expect_is(aa$setName, "character")
})

test_that("pg_list_records() works", {
  skip_on_cran()

  aa <- pg_list_records(from='2015-09-01', until='2015-09-10')

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_match(aa$identifier, "oai:pangaea.de")
  expect_is(aa$title, "character")
})

test_that("pg_list_metadata_formats() works", {
  skip_on_cran()

  aa <- pg_list_metadata_formats()

  expect_is(aa, "data.frame")
  expect_named(aa, c('metadataPrefix', 'schema', 'metadataNamespace'))
  expect_true(any(grepl("oai_dc", aa$metadataPrefix)))
})

test_that("pg_list_identifiers() works", {
  skip_on_cran()

  aa <- pg_list_identifiers(from = '2015-09-01', until = '2015-09-05')

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_match(aa$identifier, "oai:pangaea.de")
})

test_that("pg_get_record() works", {
  skip_on_cran()

  aa <- pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")

  expect_is(aa, "data.frame")
  expect_is(aa, "oai_df")
  expect_match(aa$identifier, "oai:pangaea.de")
  expect_match(aa$datestamp, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
})

test_that("fails well", {
  skip_on_cran()

  expect_error(pg_list_sets(as = "adsff"), "not in acceptable")
  expect_error(pg_list_records(prefix = "adsfadf"), "unknown")
  expect_error(pg_list_identifiers(from = 3), "Invalid datestamp")
  expect_error(pg_get_record(identifier = 4444), "Identifier is not a valid OAI")
})

##
## File:   test_05.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("XML files")

## Note the data path.
path <- system.file("extdata", "data.xml", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(2))

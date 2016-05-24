##
## File:   test_01.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("Directory trees")

## Note the data path.
path <- system.file("extdata", ".", package = "SchemaOnRead")

## Perform an overall test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(9))

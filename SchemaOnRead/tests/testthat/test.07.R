##
## File:   test_07.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("RDS files")

## Note the data path.
path <- system.file("extdata", "example.rds", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(2))

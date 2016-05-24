##
## File:   test_06.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("ARFF files")

## Note the data path.
path <- system.file("extdata", "arffexample.arff", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(6))

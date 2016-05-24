##
## File:   test_08.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("CSV files")

## Note the data path.
path <- system.file("extdata", paste("dir1", .Platform$file.sep,
  "Data.csv", sep = ""), package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(3))

##
## File:   test_10.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("XLSX files")

## Note the data path.
path <- system.file("extdata", paste("dir1", .Platform$file.sep,
  "Data1.xlsx",sep = ""), package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)$Data),
  testthat::equals(3))

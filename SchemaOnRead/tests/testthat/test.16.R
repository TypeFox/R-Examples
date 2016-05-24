##
## File:   test_16.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("REC files")

## Note the data path.
path <- system.file("extdata", paste("dir2", .Platform$file.sep,
  "example.rec", sep = ""), package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(4))

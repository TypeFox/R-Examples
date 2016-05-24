##
## File:   test_15.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("PAJ files")

## Note the data path.
path <- system.file("extdata", paste("dir2", .Platform$file.sep,
  "Example.paj", sep = ""), package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(5))

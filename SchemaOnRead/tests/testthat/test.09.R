##
## File:   test_09.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("DIF files")

## Note the data path.
path <- system.file("extdata", paste("dir1", .Platform$file.sep,
  "Data1.dif", sep = ""), package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(3))

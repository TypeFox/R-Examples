##
## File:   test_12.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("TXT files")

## Note the data path.
path <- system.file("extdata", paste("dir1", .Platform$file.sep,
  "example.txt", sep = ""), package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)$V1),
  testthat::equals(4))

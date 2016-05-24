##
## File:   test_19.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("TIFF files")

## Note the data path.
path <- system.file("extdata", "image.tiff", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(23680))

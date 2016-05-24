##
## File:   test_04.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("PNG files")

## Note the data path.
path <- system.file("extdata", "image.png", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(23680))

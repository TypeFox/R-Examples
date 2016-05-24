##
## File:   test_18.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("GIF files")

## Note the data path.
path <- system.file("extdata", "image.gif", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(4))

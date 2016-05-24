##
## File:   test_03.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("JPG files")

## Note the data path.
path <- system.file("extdata", "image.jpg", package = "SchemaOnRead")

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path)),
  testthat::equals(17760))

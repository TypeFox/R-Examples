##
## File:   test_02.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("Nonexistant entries")

## Perform a test.
testthat::expect_that(SchemaOnRead::schemaOnRead("nonexistant"),
  testthat::equals("Entry Does Not Exist"))

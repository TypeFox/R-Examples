##
## File:   test_20.R
## Author: Michael J. North
## Date:   November 25, 2015
##

## Note the type of test.
context("User processors")

## Note the data path.
path <- system.file("", "Data.xyz", package = "SchemaOnRead")

## Define a new processor.
newProcessor <- function(path, processors, verbose) {

  # Check the file existance and extensions.
  if (!SchemaOnRead::checkExtensions(path, c("xyz"))) return(NULL)

  ## As an example, attempt to read an XYZ file as a CSV file.
  read.csv(path, header = FALSE)

}

## Define a new processors list.
newProcessors <- c(newProcessor, SchemaOnRead::defaultProcessors())

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path = path,
  processors = newProcessors)), testthat::equals(1))

## Define a new processors list.
newProcessors <- c(newProcessor, SchemaOnRead::simpleProcessors())

## Perform a test.
testthat::expect_that(length(SchemaOnRead::schemaOnRead(path = path,
  processors = newProcessors)), testthat::equals(1))

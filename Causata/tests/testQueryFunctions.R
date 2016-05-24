# Test query helper functions
# Author: David Barker
###############################################################################

library(testthat)
library(Causata)

context("QueryFunctions")

source("utils.R")
equals <- testthat::equals

test_that("Quotes are escaped in sql strings", {
  expect_that(EqualTo("foo'bar'baz")$operand, equals("'foo\\'bar\\'baz'"))
})

test_that("In also escapes", {
  expect_that(In("foo'bar", "baz")$operand, equals("('foo\\'bar','baz')"))
})
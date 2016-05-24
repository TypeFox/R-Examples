# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'display-methods\'')
test_that('as.character() works', {
  res <- as.character(randTree(10))
  expect_that(class(res)[[1]], equals('character'))
})
test_that('show() works', {
  (randTree(10))
})
test_that('str() works', {
  str(randTree(10))
})
test_that('print() works', {
  print(randTree(10))
})
test_that('viz() works', {
  print(randTree(10))
})
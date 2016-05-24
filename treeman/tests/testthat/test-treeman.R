# LIBS
library(treeman)
library(testthat)

# DATA
tree <- randTree(100)

# RUNNING
context('Testing \'TreeMan Class\'')
test_that('vaildObject() works', {
  res <- validObject(tree)
  expect_that(res, is_true())
  tree@ndlst[['n2']][['id']] <- 'oh oh.... invalid ID'
  expect_error(validObject(tree))
})
test_that('[[ works', {
  nd <- sample(names(tree@ndlst), 1)
  nd <- tree[[nd]]
  expect_that(class(nd)[[1]], equals('Node'))
})
test_that('[ works', {
  expect_error(tree['not a valid slot name'])
  expect_that(tree['ntips'], equals(100))
  expect_that(tree['nnds'], equals(99))
  expect_that(tree['age'], equals(tree@age))
})
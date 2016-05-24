# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'gen-methods\'')
test_that('randTree() works', {
  ns <- sample(2:100, 5)
  for(n in ns) {
    tree <- randTree(n)
    expect_that(tree['ntips'], equals(n))
    expect_that(tree['pd'], is_more_than(tree['age']))
  }
})
# LIBS
library(treeman)
library(testthat)

# RUNNING
context('Testing \'Node Class\'')
test_that('.newNd works', {
  nd <- treeman:::.newNd(randTree(10), 'n1')
  expect_that(class(nd)[[1]], equals('Node'))
})
test_that('Node works (rooted + with spns)', {
  tree <- randTree(10)
  nid <- names(tree@ndlst)[sample(1:10, 1)]
  nd <- tree[[nid]]
})
test_that('Node works (rooted + w/o spns)', {
  tree <- randTree(10)
  tree <- setNdsSpn(tree, ids=NULL, vals=NULL)
  nid <- names(tree@ndlst)[sample(1:10, 1)]
  nd <- tree[[nid]]
})
test_that('Node works (unrooted  + w/o spns))', {
  tree <- randTree(10)
  #TODO: remove branch lengths
  #TODO: remove root
  nid <- names(tree@ndlst)[sample(1:10, 1)]
  nd <- tree[[nid]]
})
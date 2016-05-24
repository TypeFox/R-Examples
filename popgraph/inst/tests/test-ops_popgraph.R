context("ops_popgraph.R")

test_that("tests", {
  
  e1 <- as.popgraph( graph.atlas(716) )
  e2 <- as.popgraph( graph.atlas(806) )
  
  expect_that( e1 - 2, throws_error() )
  expect_that( e1 - "A", throws_error() )
  expect_that( e1 - FALSE, throws_error() )
  
  e3 <- e1 - e2
  expect_that( inherits(e3,"popgraph"), is_true() )
  expect_that( length(V(e3)), equals(7) )
  expect_that( length(E(e3)), equals(3) )
  
})
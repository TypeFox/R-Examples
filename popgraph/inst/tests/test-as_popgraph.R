context("as_popgraph.R")

test_that("testing", {
  
  A <- matrix(0,nrow=4,ncol=4)
  A[1,2] <- A[1,3] <- A[2,3] <- A[3,4] <- 1
  A <- A + t(A)
  
  graph <- as.popgraph(A)
  
  expect_that( graph, is_a("igraph") )
  expect_that( graph, is_a("popgraph"))
  expect_that( length(V(graph)), equals(4) )
  
  expect_that( as.popgraph(FALSE), throws_error() )
}
)

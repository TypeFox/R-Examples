library(igraph)
context("igraph has some broken functionality which should be fixed in the
        future but we have to rely on now. We need to make sure it works as
        expected.")

test_that("restricted argument to graph.bfs is 0 indexed for some reason", {
  g = graph.edgelist(matrix(c(1,2,2,3), ncol=2, byrow=T))
  expect_equal(as.numeric(graph.bfs(g, 1, restricted = c(1,2), unreachable=F)$order), as.numeric(c(NA,NA,NA)))
  expect_equal(as.numeric(graph.bfs(g, 1, restricted = c(0,1), unreachable=F)$order), c(1,2,NA))
})

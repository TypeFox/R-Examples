
context("Optimize node sizes")

test_that("optimize_sizes", {

  nodes_edges <- nodes_edges()

  result <- optimize_sizes(nodes_edges$nodes, nodes_edges$edges)

  expected <- structure(
    rep(1, nrow(nodes_edges$nodes)),
    names = nodes_edges$nodes[,1]
  )
  expected["draw.nodes"] <- 3
  expected["sankey"] <- 15
  expected["ypos_present"] <- 2
  expected["par"] <- 2

  expect_equal(result, expected)
})

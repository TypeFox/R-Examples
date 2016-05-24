
context("Automatic x coords")

test_that("optimize_x works", {

  nodes_edges <- nodes_edges()

  result <- optimize_x(nodes_edges$nodes, nodes_edges$edges)

  expected <- structure(
    rep(2, nrow(nodes_edges$nodes)),
    names = nodes_edges$nodes[,1]
  )
  expected["plot.sankey"] <- 0
  expected["sankey"] <- 1
  expected[c("color_ramp_palette_alpha", "ypos_present", "curveseg",
             "points", "rect", "text")] <- 3

  expect_equal(result, expected)
})

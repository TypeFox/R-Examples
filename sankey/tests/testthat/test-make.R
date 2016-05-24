
context("Making sankey plots")

test_that("make_sankey", {

  ne <- nodes_edges()
  result <- make_sankey(ne$nodes, ne$edges, y = "simple")

  expect_true(is.data.frame(result$nodes))
  expect_true(is.data.frame(result$edges))
  expect_true(all(
    c("name", "col", "size", "x", "shape", "lty", "srt", "textcol",
      "label", "adjx", "boxw", "cex", "bottom", "top", "center", "pos",
      "textx", "texty") %in% names(result$nodes)
  ))
  expect_true(all(
    c("from", "to", "colorstyle", "curvestyle",
      "col", "weight") %in% names(result$edges)
  ))

})

test_that("make_sankey handles nodes = NULL", {

  ne <- nodes_edges()
  ne$nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = sort(ne$nodes[,1])
  )

  result1 <- make_sankey(ne$nodes, ne$edges, y = "simple")
  result2 <- make_sankey(edges = ne$edges, y = "simple")

  expect_equal(result1, result2)
})

test_that("breaking edges works", {

  ne <- nodes_edges()
  nodes <- ne$nodes
  edges <- ne$edges

  nodes$col     <- nodes$col     %||% color_nodes(nodes, edges)
  nodes$shape   <- nodes$shape   %||% "rectangle"
  nodes$lty     <- nodes$lty     %||% 1
  nodes$srt     <- nodes$srt     %||% 0
  nodes$textcol <- nodes$textcol %||% "black"
  nodes$label   <- nodes$label   %||% nodes[,1]
  nodes$adjx    <- nodes$adjx    %||% 1/2
  nodes$adjy    <- nodes$adjy    %||% 1
  nodes$boxw    <- nodes$boxw    %||% 0.2
  nodes$cex     <- nodes$cex     %||% 0.7
  nodes$size    <- nodes$size    %||% optimize_sizes(nodes, edges)
  nodes$x       <- nodes$x       %||% optimize_x(nodes, edges)

  edges$colorstyle <- edges$colorstyle %||% "gradient"
  edges$curvestyle <- edges$curvestyle %||% "sin"
  edges$col        <- edges$col        %||% color_edges(nodes, edges)
  edges$weight     <- edges$weight     %||% 1

  res <- do_break_edges(nodes, edges)

  expect_equal(sort(names(res$nodes)), sort(names(nodes)))
  expect_equal(sort(names(res$edges)), sort(names(edges)))
})

test_that("breaking edges if requested", {

  ne <- nodes_edges()
  nodes <- ne$nodes
  edges <- ne$edges

  x <- make_sankey(nodes, edges, break_edges = TRUE)
  expect_equal(
    nrow(x$nodes),
    nrow(nodes) + 1
  )

})

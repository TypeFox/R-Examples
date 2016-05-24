
context("Sankey plot")

## We cannot really test this properly until we'll have
## a robust image fingerprinting package on CRAN.

test_that("sankey plots nicely", {

  ne <- nodes_edges()
  x <- make_sankey(ne$nodes, ne$edges, break_edges = FALSE)
  png(tmp <- tempfile())
  sankey(x)
  dev.off()

  expect_true(TRUE)
})

test_that("plot.sankey plots nicely", {

  ne <- nodes_edges()
  x <- make_sankey(ne$nodes, ne$edges, break_edges = FALSE)
  png(tmp <- tempfile())
  plot(x)
  dev.off()

  expect_true(TRUE)
})

test_that("sankey plots with points as nodes", {

  ne <- nodes_edges()
  ne$nodes$shape <- "point"

  x <- make_sankey(ne$nodes, ne$edges, break_edges = FALSE)
  png(tmp <- tempfile())
  plot(x)
  dev.off()

  expect_true(TRUE)
})

test_that("better optimizing edge placement", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    name = c("A", "B", "C", "D", "E", "F", "G")
  )
  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("B", "B", "B", "B", "B", "D", "D"),
    to   = c("F", "G", "C", "D", "E", "F", "G")
  )

  x <- make_sankey(nodes, edges, break_edges = TRUE, gravity = "center")
  png(tmp <- tempfile())
  plot(x)
  dev.off()

  expect_true(TRUE)
})

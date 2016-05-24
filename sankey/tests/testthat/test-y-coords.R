
context("Automatic y coords")

test_that("optimize_y", {

  ne <- nodes_edges()
  ne$nodes$x    <- optimize_x    (ne$nodes, ne$edges)
  ne$nodes$size <- optimize_sizes(ne$nodes, ne$edges)

  result <- optimize_y_simple(ne$nodes, ne$edges)

  expect_equal(result[,1:3], ne$nodes)
  expect_equal(
    result$bottom,
    c(-72.6, -66.8, -61, -55.2, -49.4, -41.6, -35.8, 0, 0, -30, -23.2,
      -17.4, -11.6, -5.8, 0, -30, -24.2, -17.4, -11.6, -5.8, 0)
  )
  expect_equal(
    result$top,
    c(-73.6, -67.8, -62, -56.2, -50.4, -44.6, -36.8, -1, -15, -31,
      -25.2, -18.4, -12.6, -6.8, -1, -31, -25.2, -19.4, -12.6, -6.8, -1)
  )
  expect_equal(
    result$center,
    c(-73.1, -67.3, -61.5, -55.7, -49.9, -43.1, -36.3, -0.5, -7.5,
      -30.5, -24.2, -17.9, -12.1, -6.3, -0.5, -30.5, -24.7, -18.4,
      -12.1, -6.3, -0.5)
  )
})

test_that("optimize_y", {

  ne <- nodes_edges()
  ne$nodes$x    <- optimize_x    (ne$nodes, ne$edges)
  ne$nodes$size <- optimize_sizes(ne$nodes, ne$edges)

  result <- optimize_y(ne$nodes, ne$edges, mode = "optimal", gravity = "top")
  result$y <- result$center
  result <- set_integer_y(result)

  xpos <- sort(unique(result$x))
  for (x in xpos) {
    mynodes <- which(result$x == x)
    for (n1 in mynodes) {
      for (n2 in mynodes) {
        if (n1 == n2) next;
        e1 <- which(ne$edges[,1] == result[n1, 1])
        e2 <- which(ne$edges[,1] == result[n2, 1])
        expect_equal(crossing_edges(result, ne$edges, e1, e2), 0)
      }
    }
  }
})

test_that("optimize_y (simple) uses supplied y coordinates", {

  ne <- nodes_edges()
  ne$nodes$x    <- optimize_x    (ne$nodes, ne$edges)
  ne$nodes$size <- optimize_sizes(ne$nodes, ne$edges)

  result <- optimize_y(
    ne$nodes,
    ne$edges,
    mode = "optimal",
    gravity = "center"
  )
  result$y <- result$center

  result <- result[ sample.int(nrow(result)), ]

  result2 <- optimize_y(
    result,
    ne$edges,
    mode = "simple",
    gravity = "top"
  )

  expect_equal(result, result2)

  result3 <- optimize_y(
    result,
    ne$edges,
    mode = "simple",
    gravity = "bottom"
  )

  expect_equal(result, result3)

  result4 <- optimize_y(
    result,
    ne$edges,
    mode = "simple",
    gravity = "center"
  )

  expect_equal(result, result4)
})

test_that("optomize_y (simple) top, bottom and center", {

  ne <- nodes_edges()
  ne$nodes$x    <- optimize_x    (ne$nodes, ne$edges)
  ne$nodes$size <- optimize_sizes(ne$nodes, ne$edges)

  top <- optimize_y_simple(ne$nodes, ne$edges, gravity = "top")
  cen <- optimize_y_simple(ne$nodes, ne$edges, gravity = "center")
  bot <- optimize_y_simple(ne$nodes, ne$edges, gravity = "bottom")

  ## TODO: how to test these?
})

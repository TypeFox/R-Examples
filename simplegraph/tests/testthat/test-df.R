
context("DF graphs")

test_that("we can create a graph from a data frame", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "c", "a", "b", "e"),
    to   = c("b", "a", "c", "e", "d", "a")
  )

  g <- graph(nodes, edges)
  expect_true(is(g, "simplegraph"))
  expect_true(is(g, "simplegraph_df"))
})

test_that("graphs without edges are OK", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = character(),
    to   = character()
  )

  g <- graph(nodes, edges)
  expect_true(is(g, "simplegraph"))
  expect_true(is(g, "simplegraph_df"))
})

test_that("graphs without nodes are OK", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = character()
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = character(),
    to   = character()
  )

  g <- graph(nodes, edges)
  expect_true(is(g, "simplegraph"))
  expect_true(is(g, "simplegraph_df"))
})

test_that("sanitize will catch problems", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "c", "a", "b", "e"),
    to   = c("b", "a", "c", "e", "d", "a")
  )

  g <- graph(nodes, edges)
  g$nodes <- NULL
  expect_error(sanitize(g), "No vertices")

  ## ---------------------

  g <- graph(nodes, edges)
  g$nodes <- data.frame(stringsAsFactors = FALSE)
  expect_error(sanitize(g), "must contain ids")

  ## ---------------------

  g <- graph(nodes, edges)
  g$nodes <- data.frame(id = 1:5)
  expect_error(sanitize(g), "First column must contain [(]character[)]")

  ## ---------------------

  g <- graph(nodes, edges)
  g$edges <- NULL
  expect_error(sanitize(g), "No edges found")

  ## ---------------------

  g <- graph(nodes, edges)
  g$edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "c")
  )
  expect_error(sanitize(g), "First two columns must contain .* ids in edge")

  ## ---------------------

  g <- graph(nodes, edges)
  g$edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b"),
    to = 1:2
  )
  expect_error(sanitize(g), "First two columns must contain .* ids in edge")

  ## ---------------------

  g <- graph(nodes, edges)
  g$edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "c", "x", "b", "e"),
    to   = c("b", "a", "c", "e", "d", "a")
  )
  expect_error(sanitize(g), "Unknown vertex id in edge data frame")
})

test_that("conversion from adj list is ok", {

  adjlist <- list(
    a = c("b", "e", "d"),
    b = c("a", "b"),
    c = character(),
    d = c("a"),
    e = c("b", "d")
  )

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    name = letters[1:5]
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "a", "a", "b", "b", "d", "e", "e"),
    to   = c("b", "e", "d", "a", "b", "a", "b", "d")
  )

  g <- graph(adjlist)

  g2 <- as_graph_data_frame(g)
  expect_equal(g2$nodes, nodes)
  expect_equal(g2$edges, edges)

  g3 <- as_graph_data_frame(g2)
  expect_equal(g2, g3)
})

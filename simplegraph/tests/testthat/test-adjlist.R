
context("Adjacency list")

test_that("we can create adjacency lists", {

  g <- graph(list(
    a = c("b", "d", "e"),
    b = c("a", "c"),
    c = character(),
    d = c("a", "b"),
    e = c("e")
  ))

  expect_true(is(g, "simplegraph_adjlist"))
  expect_true(is(g, "simplegraph"))
  expect_true(is(g, "list"))
})

test_that("we can create empty adjacency lists", {

  g <- graph(list(
    a = character(),
    b = character()
  ))

  expect_equal(order(g), 2)
  expect_equal(size(g), 0)

  g <- graph(structure(list(), names = character()))

  expect_true(is(g, "simplegraph_adjlist"))
  expect_true(is(g, "simplegraph"))
  expect_true(is(g, "list"))
})

test_that("sanitize catches stuff", {

  expect_error(
    graph(list(c("a", "b"), "a")),
    "must be named"
  )

  l <- list(a = "b", b = "a")
  names(l) <- c("", "b")
  expect_error(graph(l), "Names must be non-empty")

  l <- list(a = "b", b = "a")
  names(l) <- c("b", "b")
  expect_error(graph(l), "Duplicated names")

  expect_error(
    graph(list(a = "a", b = 1:2)),
    "must contain character vectors"
  )

  expect_error(
    graph(list(a = "b", b = c("a", "x"))),
    "Unknown vertices in adjacency list"
  )
})

test_that("conversion from data frame", {

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
  g2 <- as_graph_adjlist(g)
  expect_true(is(g2, "simplegraph_adjlist"))

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
  g2 <- as_graph_adjlist(g)
  expect_true(is(g2, "simplegraph_adjlist"))

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
  g2 <- as_graph_adjlist(g)
  expect_true(is(g2, "simplegraph_adjlist"))
})


context("Multigraphs")

test_that("is_simple", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()
  g4 <- g4()

  expect_false(is_simple(g1))
  expect_true(is_simple(g2))
  expect_true(is_simple(g3))
  expect_false(is_simple(g4))

  g1 <- as_graph_adjlist(g1)
  g2 <- as_graph_adjlist(g2)
  g3 <- as_graph_adjlist(g3)
  g4 <- as_graph_adjlist(g4)

  expect_false(is_simple(g1))
  expect_true(is_simple(g2))
  expect_true(is_simple(g3))
  expect_false(is_simple(g4))
})

test_that("is_loopy", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()
  g4 <- g4()

  expect_true(is_loopy(g1))
  expect_false(is_loopy(g2))
  expect_false(is_loopy(g3))
  expect_false(is_loopy(g4))

  g1 <- as_graph_adjlist(g1)
  g2 <- as_graph_adjlist(g2)
  g3 <- as_graph_adjlist(g3)
  g4 <- as_graph_adjlist(g4)

  expect_true(is_loopy(g1))
  expect_false(is_loopy(g2))
  expect_false(is_loopy(g3))
  expect_false(is_loopy(g4))
})

test_that("is_multigraph", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()
  g4 <- g4()

  expect_false(is_multigraph(g1))
  expect_false(is_multigraph(g2))
  expect_false(is_multigraph(g3))
  expect_true(is_multigraph(g4))

  g1 <- as_graph_adjlist(g1)
  g2 <- as_graph_adjlist(g2)
  g3 <- as_graph_adjlist(g3)
  g4 <- as_graph_adjlist(g4)

  expect_false(is_multigraph(g1))
  expect_false(is_multigraph(g2))
  expect_false(is_multigraph(g3))
  expect_true(is_multigraph(g4))
})

test_that("simplify on df", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()
  g4 <- g4()

  expect_true(is_simple(simplify(g1)))
  expect_true(is_simple(simplify(g2)))
  expect_true(is_simple(simplify(g3)))
  expect_true(is_simple(simplify(g4)))

  expect_equal(
    simplify(g1),
    graph(
      g1$nodes,
      data.frame(
        stringsAsFactors = FALSE,
        from = c("a", "b", "a", "b", "e"),
        to   = c("b", "a", "e", "d", "a")
      )
    )
  )

  expect_equal(
    simplify(g4),
    graph(
      g1$nodes,
      data.frame(
        stringsAsFactors = FALSE,
        from = c("a", "b", "a", "b", "e"),
        to   = c("b", "a", "e", "d", "a")
      )
    )
  )

})

test_that("simplify on adjlist", {

  g1 <- as_graph_adjlist(g1())
  g2 <- as_graph_adjlist(g2())
  g3 <- as_graph_adjlist(g3())
  g4 <- as_graph_adjlist(g4())

  expect_true(is_simple(simplify(g1)))
  expect_true(is_simple(simplify(g2)))
  expect_true(is_simple(simplify(g3)))
  expect_true(is_simple(simplify(g4)))

  expect_equal(
    simplify(g1),
    graph(list(
      a = c("b", "e"),
      b = c("a", "d"),
      c = character(),
      d = character(),
      e = "a"
    ))
  )

  expect_equal(
    simplify(g4),
    graph(list(
      a = c("b", "e"),
      b = c("a", "d"),
      c = character(),
      d = character(),
      e = "a"
    ))
  )
})

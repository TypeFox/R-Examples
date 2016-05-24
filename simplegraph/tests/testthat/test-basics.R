
context("Basics")

test_that("vertex_ids", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()

  expect_equal(vertex_ids(g1), letters[1:5])
  expect_equal(vertex_ids(g2), letters[1:5])
  expect_equal(vertex_ids(g3), character())
})

test_that("vertices", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()

  expect_equal(
    vertices(g1),
    data_frame(id = c("a", "b", "c", "d", "e"))
  )

  expect_equal(
    vertices(g2),
    data_frame(id = c("a", "b", "c", "d", "e"))
  )

  expect_equal(
    vertices(g3),
    data_frame(id = character())
  )
})

test_that("edges", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()

  expect_equal(edges(g1), g1$edges)
  expect_equal(edges(g2), g2$edges)
  expect_equal(edges(g3), g3$edges)
})

test_that("order and size", {

  g1 <- g1()
  g2 <- g2()
  g3 <- g3()

  expect_equal(order(g1), 5)
  expect_equal(order(g2), 5)
  expect_equal(order(g3), 0)

  expect_equal(size(g1), 6)
  expect_equal(size(g2), 0)
  expect_equal(size(g3), 0)
})

test_that("adjacent vertices", {
  g1 <- g1()
  g2 <- g2()
  g3 <- g3()

  expect_equal(
    adjacent_vertices(g1),
    list(
      a = c("b", "e", "b", "e"),
      b = c("a", "d", "a"),
      c = c("c", "c"),
      d = c("b"),
      e = c("a", "a")
    )
  )

  expect_equal(
    adjacent_vertices(g2),
    list(
      a = character(), b = character(), c = character(),
      d = character(), e = character()
    )
  )

  expect_equal(
    adjacent_vertices(g3),
    structure(list(), names = character())
  )
})

test_that("predecessors and successors", {

  G <- graph(list(A = c("B", "C"), B = "C", C = "A"))

  expect_equal(
    predecessors(G),
    list(A = c("C"), B = "A", C = c("A", "B"))
  )

  expect_equal(
    successors(G),
    list(A = c("B", "C"), B = "C", C = "A")
  )

})

test_that("incident edges", {

  G <- graph(list(A = c("B", "C"), B = "C", C = "A"))

  expect_equal(
    incident_edges(G),
    list(
      A = data_frame(from = c("A", "A"), to = c("B", "C")),
      B = data_frame(from = "B", to ="C"),
      C = data_frame(from = "C", to = "A")
    )
  )

  expect_equal(
    incident_edges(G, mode = "out"),
    list(
      A = data_frame(from = c("A", "A"), to = c("B", "C")),
      B = data_frame(from = "B", to ="C"),
      C = data_frame(from = "C", to = "A")
    )
  )

  expect_equal(
    incident_edges(G, mode = "in"),
    list(
      A = data_frame(from = "C", to = "A"),
      B = data_frame(from = "A", to = "B"),
      C = data_frame(from = c("A", "B"), to = c("C", "C"))
    )
  )

  expect_equal(
    incident_edges(G, mode = "all"),
    list(
      A = data_frame(from = c("A", "A", "C"), to = c("B", "C", "A")),
      B = data_frame(from = c("A", "B"), to = c("B", "C")),
      C = data_frame(from = c("A", "B", "C"), to = c("C", "C", "A"))
    )
  )
})

test_that("degree", {

  G <- graph(list(A = c("B", "C"), B = "C", C = "A"))

  expect_equal(
    degree(G, mode = "out"),
    c(A = 2, B = 1, C = 1)
  )

  expect_equal(
    degree(G, mode = "in"),
    c(A = 1, B = 1, C = 2)
  )

  expect_equal(
    degree(G, mode = "total"),
    c(A = 3, B = 2, C = 3)
  )

  expect_equal(
    degree(G, mode = "all"),
    c(A = 3, B = 2, C = 3)
  )
})


context("Weighted graphs")

test_that("is_weighted", {

  L <- graph(list(A = c("B", "C"), B = "C", C = "A"))
  expect_false(is_weighted(L))

  G <- graph(
    data.frame(
      stringsAsFactors = FALSE,
      id = c("a", "b", "c", "d")
    ),
    data.frame(
      stringsAsFactors = FALSE,
      from   = c("a", "a", "b", "b", "c"),
      to     = c("b", "d", "d", "c", "a"),
      weight = c( 1 ,  2 ,  1 ,  3 ,  2 )
    )
  )
  expect_true(is_weighted(G))

  G2 <- graph(
    data.frame(
      stringsAsFactors = FALSE,
      id = c("a", "b", "c", "d")
    ),
    data.frame(
      stringsAsFactors = FALSE,
      from   = c("a", "a", "b", "b", "c"),
      to     = c("b", "d", "d", "c", "a")
    )
  )
  expect_false(is_weighted(G2))

})

test_that("strength", {

  L <- graph(list(A = c("B", "C"), B = "C", C = "A"))
  expect_equal(strength(L, mode = "out"),   degree(L, mode = "out"))
  expect_equal(strength(L, mode = "in"),    degree(L, mode = "in"))
  expect_equal(strength(L, mode = "all"),   degree(L, mode = "all"))
  expect_equal(strength(L, mode = "total"), degree(L, mode = "total"))

  G <- graph(
    data.frame(
      stringsAsFactors = FALSE,
      id = c("a", "b", "c", "d")
    ),
    data.frame(
      stringsAsFactors = FALSE,
      from   = c("a", "a", "b", "b", "c"),
      to     = c("b", "d", "d", "c", "a"),
      weight = c( 1 ,  2 ,  1 ,  3 ,  2 )
    )
  )
  expect_equal(strength(G, mode = "out"),   c(a = 3, b = 4, c = 2, d = 0))
  expect_equal(strength(G, mode = "in"),    c(a = 2, b = 1, c = 3, d = 3))
  expect_equal(strength(G, mode = "all"),   c(a = 5, b = 5, c = 5, d = 3))
  expect_equal(strength(G, mode = "total"), c(a = 5, b = 5, c = 5, d = 3))

  G2 <- graph(
    data.frame(
      stringsAsFactors = FALSE,
      id = c("a", "b", "c", "d")
    ),
    data.frame(
      stringsAsFactors = FALSE,
      from   = c("a", "a", "b", "b", "c"),
      to     = c("b", "d", "d", "c", "a")
    )
  )
  strength(G2)
  expect_equal(strength(G2, mode = "out"),   degree(G2, mode = "out"))
  expect_equal(strength(G2, mode = "in"),    degree(G2, mode = "in"))
  expect_equal(strength(G2, mode = "all"),   degree(G2, mode = "all"))
  expect_equal(strength(G2, mode = "total"), degree(G2, mode = "total"))
})


context("Transpose")

test_that("transposing a df", {

  nodes <- data.frame(
    stringsAsFactors = FALSE,
    id = letters[1:5],
    col = 1:5
  )

  edges <- data.frame(
    stringsAsFactors = FALSE,
    from = c("a", "b", "c", "a", "b", "e"),
    to   = c("b", "a", "c", "e", "d", "a"),
    col  = 11:16
  )

  g <- graph(nodes, edges)
  expect_equal(g, transpose(transpose(g)))
})

test_that("transposing an adjlist", {

  g <- graph(list(
    a = c("b", "d", "e"),
    b = c("a", "c"),
    c = character(),
    d = c("a", "b"),
    e = c("e")
  ))

  expect_equal(g, transpose(transpose(g)))
})

test_that("transposing an empty adjlist", {

  g <- graph(list(
    a = character()
  ))

  expect_equal(g, transpose(transpose(g)))

  g <- graph(structure(list(), names = character()))

  expect_equal(g, transpose(transpose(g)))
})

context("convertRowsToList")

test_that("convertRowsToList", {
  expect_equal(
    convertRowsToList(matrix(1:4, 2, byrow=TRUE), as.vector = TRUE),
    list(c(1, 2), c(3, 4))
  )
  expect_equal(
    convertRowsToList(matrix(1:4, 2, byrow=TRUE), as.vector = FALSE),
    list(list(1, 2), list(3, 4))
  )
  expect_equal(
    convertRowsToList(setColNames(matrix(1:4, 2, byrow = TRUE), c("a", "b")),
      name.vector = TRUE, as.vector = FALSE),
    list(list(a=1, b=2), list(a=3, b=4))
  )
  expect_equal(
    convertRowsToList(setColNames(matrix(1:4, 2, byrow = TRUE), c("a", "b")),
      name.list = FALSE, as.vector = FALSE),
    list(list(1, 2), list(3, 4))
  )
  levs = c("a", "b")
  expect_equal(
    convertRowsToList(data.frame(a = 1:2, b = factor(c("a", "b"))),
      name.list = FALSE, factors.as.char = TRUE),
    list(list(1, "a"), list(2, "b"))
  )
  expect_equal(
    convertRowsToList(setRowNames(data.frame(a = 1:2, b = factor(c("a", "b"))), c("x", "y")),
      name.list = TRUE, name.vector = TRUE, factors.as.char = FALSE),
    list(x = list(a = 1, b = factor("a", levels = levs)), y = list(a = 2, b = factor("b", levels = levs)))
  )
})

test_that("convertColsToList", {
  expect_equal(
    convertColsToList(matrix(1:4, 2, byrow = FALSE), as.vector = TRUE),
    list(c(1, 2), c(3, 4))
  )
  expect_equal(
    convertColsToList(matrix(1:4, 2, byrow = FALSE), as.vector = FALSE),
    list(list(1, 2), list(3, 4))
  )
  expect_equal(
    convertColsToList(setRowNames(matrix(1:4, 2, byrow = FALSE), c("a", "b")),
      name.vector = TRUE, as.vector = FALSE),
    list(list(a = 1, b = 2), list(a = 3, b = 4))
  )
})


test_that("convertColsToList works with data.frame", {
  d1 = iris
  x1 = as.list(d1)
  expect_equal(convertColsToList(d1, factors.as.char = FALSE), x1)
  d2 = d1; d2$Species = as.character(d2$Species);
  x2 = as.list(d2)
  expect_equal(convertColsToList(d1, factors.as.char = TRUE), x2)
})


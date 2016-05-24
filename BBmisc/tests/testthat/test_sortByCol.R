context("sortByCol")

test_that("sortByCol", {
  d1 = setRowNames(data.frame(x = c(2, 3, 1), y = c("a", "c", "b")), c("1", "2", "3"))

  d2 = sortByCol(d1, "x")
  d3 = setRowNames(data.frame(x = c(1, 2, 3), y = c("b", "a", "c")), c(3, 1, 2))
  expect_equal(d2, d3)

  d2 = sortByCol(d1, "x", asc = FALSE)
  d3 = setRowNames(data.frame(x = c(3, 2, 1), y = c("c", "a", "b")), c(2, 1, 3))
  expect_equal(d2, d3)

  d2 = sortByCol(d1, c("x", "y"))
  d3 = setRowNames(data.frame(x = c(1, 2, 3), y = c("b", "a", "c")), c(3, 1, 2))
  expect_equal(d2, d3)

  d2 = sortByCol(d1, "y")
  d3 = setRowNames(data.frame(x = c(2, 1, 3), y = c("a", "b", "c")), c(1, 3, 2))
  expect_equal(d2, d3)

  # real tie breaker
  d1 = data.frame(x = c(2, 2, 1), y = c("a", "b", "c"))
  d2 = sortByCol(d1, c("x", "y"))
  d3 = data.frame(x = c(1, 2, 2), y = c("c", "a", "b"))
  expect_equal(d2, d3, check.attributes = FALSE)
  d2 = sortByCol(d1, c("x", "y"), asc = c(TRUE, FALSE))
  d3 = data.frame(x = c(1, 2, 2), y = c("c", "b", "a"))
  expect_equal(d2, d3, check.attributes = FALSE)

  # one col
  d1 = setRowNames(data.frame(x = c(1, 2)), c(1, 2))
  d2 = sortByCol(d1, "x", asc = FALSE)
  d3 = setRowNames(data.frame(x = c(2, 1)), c(2, 1))
  expect_equal(d2, d3)
})

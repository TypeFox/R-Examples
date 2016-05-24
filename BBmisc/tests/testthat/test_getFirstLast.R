context("getFirstLast")

test_that("getFirstLast", {
  expect_equal(getFirst(1:3), 1L)
  expect_equal(getLast(1:3), 3L)
  expect_equal(getFirst(list(iris, 1)), iris)
  expect_equal(getLast(list(iris, 1)), 1)

  expect_equal(getFirst(c(a=1, 2)), 1)
  expect_equal(names(getFirst(c(a=1, 2))), NULL)
  expect_equal(getLast(c(a=1, 2)), 2)
  expect_equal(names(getLast(c(a=1, 2))), NULL)
})


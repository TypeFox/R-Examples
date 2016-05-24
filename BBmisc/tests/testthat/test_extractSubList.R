context("extractSubList")

test_that("extractSubList", {
  xs = list(
    a = list(x = 1, y = "foo", z = matrix(1,1,1)),
    b = list(x = 2L, y = "bar", z = matrix(2,2,2))
  )
  expect_equal(extractSubList(xs, "x"), c(a = 1, b = 2))
  expect_equal(extractSubList(xs, "y"), c(a = "foo", b = "bar"))
  expect_equal(extractSubList(xs, "z"), list(a = matrix(1,1,1), b = matrix(2,2,2)))

  expect_equal(extractSubList(xs, "x", use.names = FALSE), c(1, 2))
  expect_equal(extractSubList(xs, "y", use.names = FALSE), c("foo", "bar"))
  expect_equal(extractSubList(xs, "z", use.names = FALSE), list(matrix(1,1,1), matrix(2,2,2)))

  expect_equal(extractSubList(list(), "x"), list())

  expect_equal(extractSubList(list(), "x", element.value = numeric(1)), numeric(0))
  expect_equal(extractSubList(list(), "y", element.value = character(1)), character(0))

  expect_equal(extractSubList(xs, "x", element.value = numeric(1)), c(a = 1, b = 2))
  expect_equal(extractSubList(xs, "y", element.value = character(1)), c(a = "foo", b = "bar"))

  xs = list(
    list(x = 1, y = "foo", z = matrix(1,1,1)),
    list(x = 2L, y = "bar", z = matrix(2,2,2))
  )
  expect_equal(extractSubList(xs, "y", use.names = TRUE), c("foo", "bar"))
  expect_equal(extractSubList(xs, "y", use.names = FALSE), c("foo", "bar"))

  expect_equal(
    extractSubList(list(list(a = 1:2), list(a = 3:4)), "a", simplify = "rows"),
    matrix(1:4, nrow = 2L, ncol = 2L, byrow = TRUE)
  )

  expect_equal(
    extractSubList(list(list(a = 1), list(a = 2)), "a", simplify = "rows"),
    matrix(1:2, nrow = 2L, ncol = 1)
  )
})




test_that("extractSubList works with repeated indexing", {
  xs = list(
    a = list(v = list(x = 1), w = list(y = "foo")),
    b = list(v = list(x = 2), w = list(y = "bar"))
  )
  expect_equal(extractSubList(xs, c("v", "x")), c(a = 1, b = 2))
  expect_equal(extractSubList(xs, c("w", "y")), c(a = "foo", b = "bar"))

  expect_equal(extractSubList(xs, c("v", "x"), element.value = numeric(1)), c(a = 1, b = 2))
  expect_equal(extractSubList(xs, c("w", "y"), element.value = character(1)), c(a = "foo", b = "bar"))


  expect_equal(extractSubList(xs, c("v", "x"), simplify = "rows", use.names = FALSE), matrix(c(1, 2), nrow = 2))
  expect_equal(extractSubList(xs, c("v", "x"), simplify = "cols", use.names = TRUE),
    setColNames(matrix(c(1, 2), nrow = 1), c("a", "b")))
})





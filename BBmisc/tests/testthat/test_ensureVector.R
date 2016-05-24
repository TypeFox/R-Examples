context("ensureVector")

test_that("ensureVector", {
  expect_equal(ensureVector("a", n = 2L), c("a", "a"))
  expect_equal(ensureVector("a", n = 2L, cl = "integer"),  "a")
  expect_equal(ensureVector(1, n = 1), c(1))

  expect_equal(ensureVector(c("a", "b"), n = 10L), c("a", "b"))

  expect_equal(ensureVector(iris, n = 1L), list(iris))
  expect_equal(ensureVector(iris, n = 2L, cl = "matrix"), iris)
  expect_equal(ensureVector(iris, n = 2L, cl = "data.frame"), list(iris, iris))
  expect_equal(ensureVector(iris, n = 2L), list(iris, iris))
  expect_equal(ensureVector(iris, n = 2L, names = c("a", "b")), list(a = iris, b = iris))
})


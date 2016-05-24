context("as.list")

test_that("as.list", {
  path = tempfile()
  f = fail(path)

  f$put(a = 1, b = 2, c = 3)
  expect_equal(f$as.list(), setNames(as.list(1:3), letters[1:3]))
  expect_equal(f$as.list("a"), setNames(as.list(1), letters[1]))


  # invalid keys and empty sets
  expect_equal(length(f$as.list(NULL)), 0)
  expect_error(f$as.list("xxx"))
  expect_equal(length(f$as.list(character(0L))), 0L)
})

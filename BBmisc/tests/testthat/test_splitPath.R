context("splitPath")

test_that("splitPath", {
  p = tempfile()
  res = splitPath(p)
  expect_true(is.list(res))
  expect_equal(names(res), c("drive", "path"))
  expect_true(length(res$drive) == as.integer(isWindows()))
  expect_true(is.character(res$path))
  expect_true(length(res$path) >= 1L)

  p = c("tmp", "foo", "", "bar")
  res = splitPath(collapse(p, "/"))
  expect_equal(tail(res$path, 3), p[-3])
  res = splitPath(collapse(p, "\\"))
  expect_equal(tail(res$path, 3), p[-3])
})

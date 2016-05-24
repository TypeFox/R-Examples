context("size")

test_that("size", {
  path = tempfile()
  f = fail(path)

  expect_equal(length(f$size()), 0)
  f$put(a = 1, b = 2, c = 3)
  expect_equal(length(f$size(("a"))), 1)
  expect_equal(length(f$size()), 3)
  expect_true(is.numeric(f$size()))

  expect_true(all(f$size() > f$size(unit="Mb")))
  expect_true(unname(is.na(f$size("d"))))
  expect_equal(f$size(character(0L)), setNames(numeric(0L), character(0L)))
})

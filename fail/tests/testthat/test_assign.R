context("assign")

test_that("assign", {
  path = tempfile()
  f = fail(path)

  f$put(a = 1, b = 2)
  f$assign("a")
  expect_true(exists("a"))
  expect_false(exists("b"))
  a = 3
  expect_true(a == 3)
  f$assign("a")
  expect_true(a == 1)
})

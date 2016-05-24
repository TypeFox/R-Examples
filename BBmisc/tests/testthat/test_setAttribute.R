context("setAttribute")

test_that("setAttribute", {
  x = 1:9
  x = setAttribute(x, "foo", "bar")
  x = setAttribute(x, "dim", c(3,3))
  expect_equal(attr(x, "foo"), "bar")
  expect_equal(nrow(x), 3)
  expect_equal(ncol(x), 3)
})

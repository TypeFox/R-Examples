context("namedList")

test_that("namedList", {
  expect_equal(namedList(), vector("list", 0))
  expect_equal(namedList("a"), list(a=NULL))
  expect_equal(namedList(c("a", "b")), list(a=NULL, b=NULL))
  expect_equal(namedList(c("a", "b"), 1), list(a=1, b=1))
  f = function(x) x^2
  expect_equal(namedList(c("a", "b"), f(2)), list(a=4, b=4))
})

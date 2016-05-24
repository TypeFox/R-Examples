context("setValue")

test_that("setValue", {
  xs1 = list(a=1, b=2)
  expect_equal(setValue(xs1, "b", 3), list(a=1, b=3))
  expect_equal(setValue(xs1, "b", NULL), list(a=1, b=NULL))
  expect_equal(setValue(xs1, "c", 3), list(a=1, b=2, c=3))
  expect_equal(setValue(xs1, c("a","b"), as.list(4:5)), list(a=4, b=5))
  expect_equal(setValue(xs1, c("b","c"), as.list(4:5)), list(a=1, b=4, c=5))
})


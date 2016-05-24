context("addClasses")

test_that("addClasses", {
  x = list(a=1)
  x = addClasses(x, "foo1")
  expect_equal(x, structure(list(a=1), class=c("foo1", "list")))
  x = addClasses(x, c("foo2", "foo3"))
  expect_equal(x, structure(list(a=1), class=c("foo2", "foo3", "foo1", "list")))
})

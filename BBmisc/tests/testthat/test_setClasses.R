context("setClasses")

test_that("setClasses", {
  x = list(a=1)
  expect_equal(setClasses(x, "foo"), structure(list(a=1), class="foo"))
  expect_equal(setClasses(x, c("foo1", "foo2")), structure(list(a=1), class=c("foo1", "foo2")))
})

context("def")

test_that("definitions", {
  expect_equal(format(make_def("a", "b")), "a=b")
  expect_equal(format(make_def("a", "b c")), "a=b c")
  expect_equal(format(make_def("a", "'b c'")), "a='b c'")
  expect_equal(format(make_def("a", '"b c"')), 'a="b c"')
})

test_that("error checking", {
  expect_error(make_def(LETTERS, letters), "character value")
  expect_error(make_def("a", letters), "character value")
})

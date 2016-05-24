context("filterNull")

test_that("filterNull", {
  li = list(1, 2, NULL, 3)
  expect_equal(filterNull(li), list(1, 2, 3))
  expect_equal(filterNull(list()), list())
  expect_error(filterNull(iris))
})

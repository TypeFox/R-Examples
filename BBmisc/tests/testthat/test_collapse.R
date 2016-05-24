context("collapse")

test_that("collapse", {
  expect_equal(collapse(1), "1")
  expect_equal(collapse(1:2), "1,2")
  expect_equal(collapse(c("a", "22")), "a,22")
})

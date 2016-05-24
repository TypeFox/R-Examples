context("factor")

test_that("combine", {
  x = factor(c("a", "b"))
  y = factor(c("b", "c"))
  expect_equal(cFactor(x,y), factor(c("a", "b", "b", "c")))
})
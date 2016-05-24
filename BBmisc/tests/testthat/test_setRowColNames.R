context("setRowColNames")

test_that("setRowColNames", {
  x = y = matrix(1:4, 2, 2)
  rownames(y) = c("a", "b")
  expect_equal(setRowNames(x, c("a", "b")), y)    
  colnames(y) = c("c", "d")
  expect_equal(setColNames(setRowNames(x, c("a", "b")), c("c", "d")), y)    
})

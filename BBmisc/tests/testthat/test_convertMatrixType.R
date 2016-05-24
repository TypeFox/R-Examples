context("convertMatrixType")

test_that("convertMatrixType", {
  a1 = matrix(1:4, 2L, 2L)
  a2 = matrix(as.character(1:4), 2L, 2L)
  expect_equal(convertMatrixType(a1, "character"), a2)
})

require("testthat")

test_that("Tests of consistency ratio", {
  comparisonMatrixValues = c("1","1/3","5",
                             "3","1","7",
                             "1/5","1/7","1")
  comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)
  matrix = pairwiseComparisonMatrix(comparisonMatrix)

  expect_equal(consistencyIndex(matrix), 0.03244379, tolerance = 1e-07)
  expect_equal(consistencyRatio(matrix, print.report = FALSE), 0.0623919, tolerance = 1e-07)
  expect_is(matrix, "PairwiseComparisonMatrix")
})


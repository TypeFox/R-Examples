require("testthat")

test_that("Tests of weights calculation", {
  comparisonMatrixValues = c("1","9","4",
                             "1/9","1","1/5",
                             "1/4","5","1")
  comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)
  matrix = pairwiseComparisonMatrix(comparisonMatrix)
  weights = calculateWeights(matrix)
  expect_equal(length(weights@weights), nrow(matrix@values))
})

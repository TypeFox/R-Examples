require("testthat")

test_that("Tests of AHP calculation", {
  comparisonMatrixValues = c("1","9","4",
                             "1/9","1","1/5",
                             "1/4","5","1")
  comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)
  matrix = pairwiseComparisonMatrix(comparisonMatrix)
  weights = calculateWeights(matrix)

  data = c(1,2,3,4,5,6,7,8,9,4,5,6)
  data = matrix(data, nrow=length(data)/3, ncol=3, byrow=TRUE)

  ahpResult = calculateAHP(weights, data)

  expect_true(ahpResult[1]<ahpResult[2])
  expect_true(ahpResult[2]<ahpResult[3])
  expect_true(ahpResult[2]==ahpResult[4])

  comparisonMatrixValues = c(1,9,4,
                             1/9,1,1/5,
                             1/4,5,1)
  comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)
  matrix = pairwiseComparisonMatrix(comparisonMatrix)

})

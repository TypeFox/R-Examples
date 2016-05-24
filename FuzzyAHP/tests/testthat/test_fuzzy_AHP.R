require("testthat")

test_that("Tests of fuzzy AHP calculation", {

  folder = "./"
  # folder = "./tests/testthat/"
  file = "testing_pairwise_comparison_matrix.txt"

  comparisonMatrix = read.csv(paste(folder, file, sep = ""), sep = ";",
                              stringsAsFactors = FALSE, header = FALSE, strip.white = TRUE)

  comparisonMatrix = as.matrix(comparisonMatrix)

  comparisonMatrix = pairwiseComparisonMatrix(comparisonMatrix)

  CR = consistencyRatio(comparisonMatrix, print.report = FALSE)

  fuzzyComparisonMatrix = fuzzyPairwiseComparisonMatrix(comparisonMatrix)

  CRF = consistencyRatio(fuzzyComparisonMatrix, print.report = FALSE)

  expect_equal(CR, CRF, tolerance = 1e-07)

  weights = calculateWeights(fuzzyComparisonMatrix)

  expect_is(fuzzyComparisonMatrix, "FuzzyPairwiseComparisonMatrix")
  expect_is(weights, "FuzzyWeights")
})


test_that("Tests of fuzzy pairwise comparsion matrix consturction from character matrix", {

  matrix = c("1;1;1", "8;9;9", "5;6;7",
              "1/9;1/9;1/8", "1;1;1", "3;4;5",
             "1/7;1/6;1/5", "1/5;1/4;1/3", "1;1;1")
  matrix = matrix(matrix, nrow=3, ncol=3, byrow = TRUE)
  fuzzyMatrix = fuzzyPairwiseComparisonMatrix(matrix)

  expect_is(fuzzyMatrix, "FuzzyPairwiseComparisonMatrix")

# ------------------------------------------------------------------
  matrix = c("1;1;1", "8;9;9",
             "1/9;1/9;1/8", "1;1;1",
             "1/7;1/6;1/5", "1/5;1/4;1/3")
  matrix = matrix(matrix, nrow=3, ncol=2, byrow = TRUE)

  expect_error(fuzzyPairwiseComparisonMatrix(matrix), "The fuzzy pairwise comparison matrix is not a square matrix.")

# ------------------------------------------------------------------
  matrix = c("1;1;1", "8;9;9", "5;6;7",
             "1/9;1/9;1/8", "1;1;1", "3;4;5",
             "1/7;1/6;1/5", "1/5;1/4;1/3", "1;1;3")
  matrix = matrix(matrix, nrow=3, ncol=3, byrow = TRUE)

  expect_error(fuzzyPairwiseComparisonMatrix(matrix), "(1,1,1)")
  }
)

## ----setup, include=FALSE---------------------------------------------------------------------------------------------
options(width=120)
library(FuzzyAHP)

## ---- eval = FALSE----------------------------------------------------------------------------------------------------
#  matrixFile =  "comparison_matrix.csv"
#  comparisonMatrix = read.csv(matrixFile, sep = ";",
#                     stringsAsFactors = FALSE, header = TRUE, row.names = 1, strip.white = TRUE)
#  comparisonMatrix = as.matrix(comparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
comparisonMatrixValues = c(1,9,5,
                           1/9,1,1/3,
                           1/5,3,1)
comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)

## ---------------------------------------------------------------------------------------------------------------------
comparisonMatrixValues = c("1","9","5",
                           "1/9","1","1/3",
                           "1/5","3","1")
comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)

## ---------------------------------------------------------------------------------------------------------------------
comparisonMatrix = pairwiseComparisonMatrix(comparisonMatrix)
print(comparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
CR = consistencyRatio(comparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
weakConsistency = weakConsistency(comparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
strictConsistency = strictConsistency(comparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
weights = calculateWeights(comparisonMatrix)
print(weights)

## ---------------------------------------------------------------------------------------------------------------------
values = c(4,5,3,
1,3,9,
8,6,4,
3,2,7,
6,7,5,
4,5,3,
NA,9,9,
NA,NA,NA)
values = matrix(values, nrow = length(values)/length(weights@weights), ncol = length(weights@weights), byrow = TRUE)

## ---------------------------------------------------------------------------------------------------------------------
result = calculateAHP(weights, values)
print(result)

## ---------------------------------------------------------------------------------------------------------------------
rank = compareResults(result)
print(rank)

## ---------------------------------------------------------------------------------------------------------------------
result = cbind(values, result, rank)
colnames(result) = c("crit1", "crit2", "crit3", "result_value", "ranking")
print(result)

## ---------------------------------------------------------------------------------------------------------------------
comparisonMatrixValues = c("1","9","5",
                       "1/9","1","1/3",
                       "1/5","3","1")
comparisonMatrix = matrix(comparisonMatrixValues, nrow = 3, ncol = 3, byrow = TRUE)
comparisonMatrix = pairwiseComparisonMatrix(comparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
fuzzyComparisonMatrix = fuzzyPairwiseComparisonMatrix(comparisonMatrix)
print(fuzzyComparisonMatrix)

## ---------------------------------------------------------------------------------------------------------------------
result = calculateAHP(fuzzyComparisonMatrix, values)

## ---------------------------------------------------------------------------------------------------------------------
fuzzyNumer = getFuzzyNumber(result, as.integer(2))
print(fuzzyNumer)

## ---------------------------------------------------------------------------------------------------------------------
defuzzified = defuzziffy(result, "Yager")
print(defuzzified)
rank = (nrow(values)+1) - rank(defuzzified, na.last = FALSE, ties.method= "max")
print(rank)

## ---------------------------------------------------------------------------------------------------------------------
ranked = compareFuzzyNumbers(result, "Chen")
print(ranked)

## ---- results = "hide"------------------------------------------------------------------------------------------------
ranked = compareFuzzyNumbers(result, "possibilityTheory")
# ranked = compareFuzzyNumbers(result, "possibilityTheory", progressBar = TRUE)

## ---------------------------------------------------------------------------------------------------------------------
print(ranked)


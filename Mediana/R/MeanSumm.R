############################################################################################################################

# Function: MeanSumm.
# Argument: Descriptive statistics across multiple simulation runs (vector or matrix), method parameters (not used in this function).
# Description: Compute mean for the vector of statistics or in each column of the matrix.

MeanSumm = function(test.result, statistic.result, parameter) {

  result = apply(statistic.result, 2, mean)
  return(result)
}
# End of MeanSumm
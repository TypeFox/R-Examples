############################################################################################################################

# Function: MedianSumm.
# Argument: Descriptive statistics across multiple simulation runs (vector or matrix), method parameters (not used in this function).
# Description: Compute median for the vector of statistics or in each column of the matrix.

MedianSumm = function(test.result, statistic.result, parameter) {

  result = apply(statistic.result, 2, stats::median)
  return(result)
}
# End of MedianSumm
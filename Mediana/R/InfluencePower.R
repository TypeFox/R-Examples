############################################################################################################################

# Function: InfluencePower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results,
#       criterion parameter (Influence cutoff).
# Description: Compute probability that the influence condition is met.

InfluencePower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: InfluencePower: alpha parameter must be specified.")
  if (is.null(parameter$cutoff)) stop("Evaluation model: InfluencePower: cutoff parameter must be specified.")

  alpha = parameter$alpha
  cutoff_influence  = parameter$cutoff

  significant = ((test.result[,1] <= alpha & test.result[,2] <= alpha & statistic.result[,1] >= cutoff_influence) | (test.result[,1] <= alpha & test.result[,2] > alpha))


  power = mean(significant)
  return(power)
}
# End of InfluencePower
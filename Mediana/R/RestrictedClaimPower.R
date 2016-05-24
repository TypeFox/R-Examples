############################################################################################################################

# Function: RestrictedClaimPower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results,
#       criterion parameter (Type I error rate and Influence cutoff).
# Description: Compute probability of restricted claim (new treatment is effective in the target population only)

RestrictedClaimPower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: RestrictedClaimPower: alpha parameter must be specified.")
  if (is.null(parameter$cutoff_influence)) stop("Evaluation model: RestrictedClaimPower: cutoff parameter must be specified.")

  alpha = parameter$alpha
  cutoff_influence  = parameter$cutoff_influence

  significant = ((test.result[,1] > alpha & test.result[,2] <= alpha) | (test.result[,1] <= alpha & test.result[,2] <= alpha & statistic.result[,1] < cutoff_influence))

  power = mean(significant)
  return(power)
}
# End of RestrictedClaimPower

############################################################################################################################

# Function: BroadClaimPower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results,
#       criterion parameter (Type I error rate and Influence cutoff).
# Description: Compute probability of broad claim (new treatment is effective in the overall population without substantial effect in the subgroup of interest)

BroadClaimPower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: BroadClaimPower: alpha parameter must be specified.")
  if (is.null(parameter$cutoff_influence)) stop("Evaluation model: BroadClaimPower: cutoff_influence parameter must be specified.")
  if (is.null(parameter$cutoff_interaction)) stop("Evaluation model: BroadClaimPower: cutoff_interaction parameter must be specified.")

  alpha = parameter$alpha
  cutoff_influence  = parameter$cutoff_influence
  cutoff_interaction  = parameter$cutoff_interaction

  significant = ((test.result[,1] <= alpha & test.result[,2] <= alpha & statistic.result[,1] >= cutoff_influence & statistic.result[,2] < cutoff_interaction) | (test.result[,1] <= alpha & test.result[,2] > alpha))

  power = mean(significant)
  return(power)
}
# End of BroadClaimPower

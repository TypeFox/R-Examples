############################################################################################################################

# Function: InteractionPower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results,
#       criterion parameter (Influence cutoff).
# Description: Compute probability that the interaction condition is met.

InteractionPower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: InteractionPower: alpha parameter must be specified.")
  if (is.null(parameter$cutoff_influence)) stop("Evaluation model: InteractionPower: cutoff_influence parameter must be specified.")
  if (is.null(parameter$cutoff_interaction)) stop("Evaluation model: InteractionPower: cutoff_interaction parameter must be specified.")

  alpha = parameter$alpha
  cutoff_influence  = parameter$cutoff_influence
  cutoff_interaction  = parameter$cutoff_interaction

  significant = (test.result[,1] <= alpha & test.result[,2] <= alpha & statistic.result[,1] >= cutoff_influence & statistic.result[,2] >= cutoff_interaction)

  power = mean(significant)
  return(power)
}
# End of InteractionPower
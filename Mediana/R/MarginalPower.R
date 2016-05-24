############################################################################################################################

# Function: MarginalPower.
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results (not used in this function),
#       criterion parameter (Type I error rate).
# Description: Compute marginal power for the vector of test results (vector of p-values or each column of the p-value matrix).

MarginalPower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: MarginalPower: alpha parameter must be specified.")

  alpha = parameter$alpha
  significant = (test.result <= alpha)
  if (is.numeric(test.result)) power = mean(significant, na.rm = TRUE)
  if (is.matrix(test.result)) power = colMeans(significant, na.rm = TRUE)

  return(power)

}
# End of MarginalPower
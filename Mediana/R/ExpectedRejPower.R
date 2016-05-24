############################################################################################################################

# Function: ExpectedRejPower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results (not used in this function),
#       criterion parameter (Type I error rate and weigth).
# Description: Compute expected number of rejected hypothesis for the test results (vector of p-values or each column of the p-value matrix).

ExpectedRejPower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: WeightedPower: alpha parameter must be specified.")

  # Get the parameter
  alpha = parameter$alpha
  ntests = ncol(test.result)
  weight = rep(1/ntests,ntests)

  significant = (test.result <= alpha)
  if (is.numeric(test.result))
    # Only one test is specified and no weight is applied
    power = mean(significant, na.rm = TRUE)
  if (is.matrix(test.result)) {
    # Weights are applied when two or more tests are specified
    # Check if the number of tests equals the number of weights
    marginal.power = colMeans(significant)
    power = ntests * sum(marginal.power * weight, na.rm = TRUE)
  }
  return(power)
}


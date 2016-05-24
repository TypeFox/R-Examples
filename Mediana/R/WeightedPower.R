############################################################################################################################

# Function: WeightedPower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results (not used in this function),
#       criterion parameter (Type I error rate and weigth).
# Description: Compute weighted power for the test results (vector of p-values or each column of the p-value matrix).

WeightedPower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: WeightedPower: alpha parameter must be specified.")
  if (is.null(parameter$weight)) stop("Evaluation model: WeightedPower: weight parameter must be specified.")
  if (length(parameter$weight) != ncol(test.result)) stop("Evaluation model: WeightedPower: The number of test weights must be equal to the number of tests.")
  if (sum(parameter$weight) != 1) stop("Evaluation model: WeightedPower: sum of weights must be equal to 1.")

  # Get the parameter
  alpha = parameter$alpha
  weight = parameter$weight

  significant = (test.result <= alpha)
  if (is.numeric(test.result))
    # Only one test is specified and no weight is applied
    power = mean(significant, na.rm = TRUE)
  if (is.matrix(test.result)) {
    # Weights are applied when two or more tests are specified
    # Check if the number of tests equals the number of weights
    marginal.power = colMeans(significant)

    power = sum(marginal.power * weight, na.rm = TRUE)
  }
  return(power)
}

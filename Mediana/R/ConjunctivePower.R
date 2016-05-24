############################################################################################################################

# Function: ConjunctivePower
# Argument: Test results (p-values) across multiple simulation runs (vector or matrix), statistic results (not used in this function),
#       criterion parameter (Type I error rate).
# Description: Compute conjunctive power for the test results (vector of p-values or each column of the p-value matrix).

ConjunctivePower = function(test.result, statistic.result, parameter) {

  # Error check
  if (is.null(parameter$alpha)) stop("Evaluation model: ConjunctivePower: alpha parameter must be specified.")

  alpha = parameter$alpha

  if (is.numeric(test.result))
    significant = (test.result <= alpha)
  if (is.matrix(test.result))
    significant = (rowSums(test.result <= alpha) == ncol(test.result))

  power = mean(significant, na.rm = TRUE)
  return(power)
}
# End of ConjunctivePower
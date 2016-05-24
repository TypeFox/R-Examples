######################################################################################################################

# Function: PropTest .
# Argument: Data set and parameter (call type and Yates' correction).
# Description: Computes one-sided p-value based on two-sample proportion test.

PropTest = function(sample.list, parameter) {
  # Determine the function call, either to generate the p-value or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Yates' correction is set up by default to FALSE
    if (is.na(parameter[[2]])) yates = FALSE
    else yates = parameter[[2]]$yates

    # Sample list is assumed to include two data frames that represent two analysis samples

    # Outcomes in Sample 1
    outcome1 = sample.list[[1]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome1.complete = outcome1[stats::complete.cases(outcome1)]

    # Outcomes in Sample 2
    outcome2 = sample.list[[2]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome2.complete = outcome2[stats::complete.cases(outcome2)]

    # One-sided p-value (treatment effect in sample 2 is expected to be greater than in sample 1)
    result = stats::prop.test(c(sum(outcome2.complete), sum(outcome1.complete)),
                              n = c(length(outcome2.complete), length(outcome1.complete)), alternative = "greater", correct = yates)$p.value
  }
  else if (call == TRUE) {

    result=list("Test for proportions")
  }

  return(result)
}
# End of PropTest
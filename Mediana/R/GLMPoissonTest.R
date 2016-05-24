######################################################################################################################

# Function: GLMPoissonTest .
# Argument: Data set and parameter (call type).
# Description: Computes one-sided p-value based on Poisson regression.
#' @importFrom stats poisson
GLMPoissonTest = function(sample.list, parameter) {

  # Determine the function call, either to generate the p-value or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Sample list is assumed to include two data frames that represent two analysis samples

    # Outcomes in Sample 1
    outcome1 = sample.list[[1]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome1.complete = outcome1[stats::complete.cases(outcome1)]

    # Outcomes in Sample 2
    outcome2 = sample.list[[2]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome2.complete = outcome2[stats::complete.cases(outcome2)]

    # Data frame
    data.complete = data.frame(rbind(cbind(2, outcome2.complete), cbind(1, outcome1.complete)))
    colnames(data.complete) = c("TRT", "RESPONSE")
    data.complete$TRT=as.factor(data.complete$TRT)

    # One-sided p-value
    result = summary(stats::glm(RESPONSE ~ TRT, data = data.complete, family = poisson))$coefficients["TRT2", "Pr(>|z|)"]/2
  }
  else if (call == TRUE) {

    result=list("Poisson regression test")
  }

  return(result)
}
# End of GLMPoissonTest
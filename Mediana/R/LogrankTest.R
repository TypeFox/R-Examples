######################################################################################################################

# Function: LogrankTest .
# Argument: Data set and parameter.
# Description: Computes one-sided p-value based on log-rank test.

LogrankTest = function(sample.list, parameter) {

  # Determine the function call, either to generate the p-value or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Sample list is assumed to include two data frames that represent two analysis samples

    # Outcomes in Sample 1
    outcome1 = sample.list[[1]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome1.complete = outcome1[stats::complete.cases(outcome1)]
    # Observed events in Sample 1 (negation of censoring indicators)
    event1 = !sample.list[[1]][, "patient.censor.indicator"]
    event1.complete = event1[stats::complete.cases(outcome1)]
    # Sample size in Sample 1
    n1 = length(outcome1.complete)

    # Outcomes in Sample 2
    outcome2 = sample.list[[2]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome2.complete = outcome2[stats::complete.cases(outcome2)]
    # Observed events in Sample 2 (negation of censoring indicators)
    event2 = !sample.list[[2]][, "patient.censor.indicator"]
    event2.complete = event2[stats::complete.cases(outcome2)]
    # Sample size in Sample 2
    n2 = length(outcome2.complete)

    # Create combined samples of outcomes, censoring indicators (all events are observed) and treatment indicators
    outcome = c(outcome1.complete, outcome2.complete)
    event = c(event1.complete, event2.complete)
    treatment = c(rep(0, n1), rep(1, n2))

    # Apply log-rank test
    surv.test = survival::survdiff(survival::Surv(outcome, event) ~ treatment)

    # Compute one-sided p-value
    result = stats::pchisq(surv.test$chisq, df = 1, lower.tail = FALSE)/2

  }

  else if (call == TRUE) {
    result = list("Log-rank test")
  }

  return(result)
}
# End of LogrankTest
######################################################################################################################

# Compute the ratio of effect sizes for HR (time-to-event) based on non-missing values in the combined sample

RatioEffectSizeEventStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=4)
      stop("Analysis model: Four samples must be specified in the RatioEffectSizeEventStat statistic.")

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

    # Get the HR from the Cox-test
    result1 = log(1 / summary(survival::coxph(survival::Surv(outcome, event) ~ treatment))$coef[,"exp(coef)"])

    # Outcomes in Sample 3
    outcome3 = sample.list[[3]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome3.complete = outcome3[stats::complete.cases(outcome3)]
    # Observed events in Sample 3 (negation of censoring indicators)
    event3 = !sample.list[[3]][, "patient.censor.indicator"]
    event3.complete = event3[stats::complete.cases(outcome3)]
    # Sample size in Sample 3
    n3 = length(outcome3.complete)

    # Outcomes in Sample 4
    outcome4 = sample.list[[4]][, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome4.complete = outcome4[stats::complete.cases(outcome4)]
    # Observed events in Sample 4 (negation of censoring indicators)
    event4 = !sample.list[[4]][, "patient.censor.indicator"]
    event4.complete = event4[stats::complete.cases(outcome4)]
    # Sample size in Sample 4
    n4 = length(outcome4.complete)

    # Create combined samples of outcomes, censoring indicators (all events are observed) and treatment indicators
    outcome = c(outcome3.complete, outcome4.complete)
    event = c(event3.complete, event4.complete)
    treatment = c(rep(0, n3), rep(1, n4))

    # Get the HR from the Cox-test
    result2 = log(1 / summary(survival::coxph(survival::Surv(outcome, event) ~ treatment))$coef[,"exp(coef)"])

    # Caculate the ratio of effect size
    result = result1 / result2

  }

  else if (call == TRUE) {
    result = list("Ratio of effect size")
  }

  return(result)
}
# End of RatioEffectSizeEventStat
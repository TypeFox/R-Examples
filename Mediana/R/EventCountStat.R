######################################################################################################################

# Compute the number of events based on non-missing values in the combined sample

EventCountStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list) == 0)
      stop("Analysis model: One sample must be specified in the EventCountStat statistic.")

    # Merge the samples in the sample list
    sample1 = do.call(rbind, sample.list)

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome1 = sample1[, "outcome"]
    # Remove the missing values due to dropouts/incomplete observations
    outcome1.complete = outcome1[stats::complete.cases(outcome1)]
    # Observed events in Sample 1 (negation of censoring indicators)
    event1 = !sample1[, "patient.censor.indicator"]
    event1.complete = event1[stats::complete.cases(outcome1)]
    # Number of events in Sample 1
    result = sum(event1.complete)

  }

  else if (call == TRUE) {
    result = list("Number of Events")
  }

  return(result)
}
# End of EventCountStat
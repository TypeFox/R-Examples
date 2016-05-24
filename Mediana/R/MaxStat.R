######################################################################################################################

# Compute the min based on non-missing values in the combined sample

MaxStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=1)
      stop("Analysis model : Only one sample must be specified in the MaxStat statistic.")

    sample = sample.list[[1]]

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome = sample[, "outcome"]
    result = max(stats::na.omit(outcome))

  }

  else if (call == TRUE) {
    result = list("Maximum")
  }

  return(result)
}
# End of MaxStat
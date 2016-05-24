######################################################################################################################

# Compute the median based on non-missing values in the combined sample

MedianStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=1)
      stop("Analysis model: Only one sample must be specified in the MedianStat statistic.")

    sample = sample.list[[1]]

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome = sample[, "outcome"]
    result = stats::median(stats::na.omit(outcome))

  }

  else if (call == TRUE) {
    result = list("Median")
  }

  return(result)
}
# End of MedianStat
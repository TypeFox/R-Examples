######################################################################################################################

# Compute the mean based on non-missing values in the combined sample

MeanStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=1)
      stop("Analysis model : Only one sample must be specified in the MeanStat statistic.")

    sample = sample.list[[1]]

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome = sample[, "outcome"]
    result = mean(stats::na.omit(outcome))

  }

  else if (call == TRUE) {
    result = list("Mean")
  }

  return(result)
}
# End of MeanStat
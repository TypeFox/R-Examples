######################################################################################################################

# Compute the difference of means between two samples for continuous variable based on non-missing values in the combined sample

DiffMeanStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=2)
      stop("Analysis model: Two samples must be specified in the DiffMeanStat statistic.")

    # Merge the samples in the sample list
    sample1 = sample.list[[1]]

    # Merge the samples in the sample list
    sample2 = sample.list[[2]]

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome1 = sample1[, "outcome"]
    outcome2 = sample2[, "outcome"]
    mean1 = mean(stats::na.omit(outcome1))
    mean2 = mean(stats::na.omit(outcome2))
    result = (mean1 - mean2)

  }

  else if (call == TRUE) {
    result = list("Difference of means")
  }

  return(result)
}
# End of DiffMeanStat
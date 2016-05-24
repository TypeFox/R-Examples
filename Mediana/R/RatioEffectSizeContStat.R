######################################################################################################################

# Compute the ratio of effect sizes for continuous based on non-missing values in the combined sample

RatioEffectSizeContStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=4)
      stop("Analysis model: Four samples must be specified in the RatioEffectSizeContStat statistic.")

    # Merge the samples in the sample list
    sample1 = sample.list[[1]]

    # Merge the samples in the sample list
    sample2 = sample.list[[2]]

    # Merge the samples in the sample list
    sample3 = sample.list[[3]]

    # Merge the samples in the sample list
    sample4 = sample.list[[4]]

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome1 = sample1[, "outcome"]
    outcome2 = sample2[, "outcome"]
    mean1 = mean(stats::na.omit(outcome1))
    mean2 = mean(stats::na.omit(outcome2))
    sdcom1 = stats::sd(c(stats::na.omit(outcome1),stats::na.omit(outcome2)))
    result1 = (mean2 - mean1) / sdcom1

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome3 = sample3[, "outcome"]
    outcome4 = sample4[, "outcome"]
    mean3 = mean(stats::na.omit(outcome3))
    mean4 = mean(stats::na.omit(outcome4))
    sdcom2 = stats::sd(c(stats::na.omit(outcome3),stats::na.omit(outcome4)))
    result2 = (mean4 - mean3) / sdcom2

    # Caculate the ratio of effect size
    result = result1 / result2

  }

  else if (call == TRUE) {
    result = list("Ratio of effect size")
  }

  return(result)
}
# End of RatioEffectSizeContStat
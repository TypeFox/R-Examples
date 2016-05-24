# Compute the ratio of effect sizes for proportions based on non-missing values in the combined sample

RatioEffectSizePropStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=4)
      stop("Analysis model: Four samples must be specified in the RatioEffectSizePropStat statistic.")

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
    prop1 = mean(stats::na.omit(outcome1))
    prop2 = mean(stats::na.omit(outcome2))
    prop = (prop2 + prop1) / 2
    result1 = (prop2 - prop1) / sqrt(prop * (1-prop))

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome3 = sample3[, "outcome"]
    outcome4 = sample4[, "outcome"]
    prop3 = mean(stats::na.omit(outcome3))
    prop4 = mean(stats::na.omit(outcome4))
    prop = (prop3 + prop4) / 2
    result2 = (prop4 - prop3) / sqrt(prop * (1-prop))

    # Caculate the ratio of effect size
    result = result1 / result2

  }

  else if (call == TRUE) {
    result = list("Ratio of effect size")
  }

  return(result)
}
# End of RatioEffectSizePropStat

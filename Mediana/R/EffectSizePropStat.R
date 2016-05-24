######################################################################################################################

# Compute the effect size for binary data based on non-missing values in the combined sample

EffectSizePropStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)!=2)
      stop("Analysis model: Two samples must be specified in the EffectSizePropStat statistic.")

    # Merge the samples in the sample list
    sample1 = sample.list[[1]]

    # Merge the samples in the sample list
    sample2 = sample.list[[2]]

    # Select the outcome column and remove the missing values due to dropouts/incomplete observations
    outcome1 = sample1[, "outcome"]
    outcome2 = sample2[, "outcome"]
    prop1 = mean(stats::na.omit(outcome1))
    prop2 = mean(stats::na.omit(outcome2))
    prop = (prop2 + prop1) / 2
    result = (prop2 - prop1) / sqrt(prop * (1-prop))

  }

  else if (call == TRUE) {
    result = list("Effect size")
  }

  return(result)
}
# End of EffectSizePropStat
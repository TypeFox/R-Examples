######################################################################################################################

# Compute the number of patients generated based on non-missing values in the combined sample

PatientCountStat = function(sample.list, parameter) {

  # Determine the function call, either to generate the statistic or to return description
  call = (parameter[[1]] == "Description")

  if (call == FALSE | is.na(call)) {

    # Error checks
    if (length(sample.list)==0)
      stop("Analysis model: One sample must be specified in the PatientCountStat statistic.")

    # Merge the samples in the sample list
    sample1 = do.call(rbind, sample.list)

    result = nrow(sample1)

  }

  else if (call == TRUE) {
    result = list("Number of Patients")
  }

  return(result)
}
# End of PatientCountStat
#######################################################################################################################

# Function: CreateDataScenarioSampleSize .
# Argument: Data frame of patients and sample size.
# Description: Create data stack for the current sample size. This function is used in the CreateDataStack function when the user uses the SampleSize.

CreateDataScenarioSampleSize = function(current.design.outcome.variables, current.sample.size) {

  # List of current data scenario
  current.data.scenario = list()
  current.data.scenario.index = 0

  # Get the number of samples
  n.samples = length(current.design.outcome.variables)

  # Get the number of outcome
  n.outcomes = length(current.design.outcome.variables[[1]])

  # For each sample, generate the data for each outcome for the current sample size
  for (sample.index in 1:n.samples){

    for (outcome.index in 1:n.outcomes){

      # Increment the index
      current.data.scenario.index = current.data.scenario.index + 1

      # Get the current sample.size for the current sample
      current.sample.size.sample = as.numeric(current.sample.size[sample.index])

      # Get the data for the current sample.size
      current.data = current.design.outcome.variables[[sample.index]][[outcome.index]]$data[(1:current.sample.size.sample),]

      # Get the sample id
      current.id = current.design.outcome.variables[[sample.index]][[outcome.index]]$id

      # Get the outcome type
      current.outcome.type = current.design.outcome.variables[[sample.index]][[outcome.index]]$outcome.type

      # Add the current sample in the list
      current.data.scenario[[current.data.scenario.index]] = list(id = current.id,
                                                                  outcome.type = current.outcome.type,
                                                                  data = current.data
      )
    }

  }

  # Return the object
  return(current.data.scenario)

} # End of CreateDataScenarioSampleSize

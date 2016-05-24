#######################################################################################################################

# Function: CreateDataScenarioEvent.
# Argument: Data frame of patients and number of events.
# Description: Create data stack for the current number of events. This function is used in the CreateDataStack function when the user uses the Event.

CreateDataScenarioEvent = function(current.design.outcome.variables, current.events, rando.ratio) {

  # List of current data scenario
  current.data.scenario = list()
  current.data.scenario.index = 0

  # Get the number of samples
  n.samples = length(current.design.outcome.variables)

  # Get the number of outcome
  n.outcomes = length(current.design.outcome.variables[[1]])

  # Get the patient indicator censor of the primary outcome from all samples
  current.design.outcome.variables.primary = lapply(current.design.outcome.variables, function(x) !x[[1]]$data[,"patient.censor.indicator"])

  # Add rows in case of unbalance randomization to bind by column
  maxrow = max(unlist(lapply(current.design.outcome.variables.primary, length)))
  current.design.outcome.variables.primary.complete = mapply(cbind,lapply(current.design.outcome.variables.primary, function(x) 	{ length(x) = maxrow
  return(x)
  }))

  # Calculate the cumulative number of events for each sample according to the randomization ratio
  n.events.cum = mapply(function(x,y) cumsum(x)[seq(y,length(x), y)], current.design.outcome.variables.primary, as.list(rando.ratio))
  index.patient = which(rowSums(n.events.cum)>=current.events)[1]

  # Get the number of patients required to get the current number of events in each sample
  index.patient = rando.ratio*index.patient

  # For each sample, generate the data for each outcome for the current sample size
  for (sample.index in 1:n.samples){

    for (outcome.index in 1:n.outcomes){

      # Increment the index
      current.data.scenario.index = current.data.scenario.index + 1

      # Get the data for the current sample.size
      current.data = current.design.outcome.variables[[sample.index]][[outcome.index]]$data[(1:index.patient[[sample.index]]),]

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

} # End of CreateDataScenarioEvent

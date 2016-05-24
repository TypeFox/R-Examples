#######################################################################################################################

# Function: CreateDataSlice.
# Argument: Data scenario (results of a single simulation run for a single combination of data scenario factores),
# list of analysis samples for defining a data slice, slice criterion and slice value.
# Description: Creates a subset of the original data set (known as a slice) based on the specified parameter. This function is useful for
# implementing interim looks and requires that the enrollment parameters should be specified in the data model.
# If parameter = "sample.size", the time point for defining the slice is determined by computing the time when X patients
# in the combined samples from the sample list complete the trial (patients who dropped out of the trial or who are censored
# are not counted as completers).
# If parameter = "event", the time point for defining the slice is determined by computing the time when X events occurs
# in the combined samples specified in the sample list.
# If parameter = "time", the data slice includes the patients who complete the trial by the specified time cutoff
# (patients who dropped out of the trial or who are censored are not counted as completers).
# After the time cutoff is determined, the data slice is then created by keeping only the
# patients who complete the trial prior to the time cutoff. The patients not included in the data slice are
# assumed to have dropped out of the trial (if the outcome type if "standard") or have been censored (if the outcome type if "event").

# X is specified by the value argument.

CreateDataSlice = function(data.scenario, sample.list, parameter, value) {

  # Number of analysis samples within the data scenario
  n.analysis.samples = length(data.scenario$sample)

  # Create the data slice
  data.slice = data.scenario

  # Select the samples from the sample list within the current data scenario
  selected.samples = list()

  index = 1
  for (i in 1:n.analysis.samples) {
    if (data.scenario$sample[[i]]$id %in% sample.list) {
      selected.samples[[index]] = data.scenario$sample[[i]]$data
      index = index + 1
    }
  }

  # Merge the selected analysis samples
  selected.analysis.sample = do.call(rbind, selected.samples)

  # Determine the time cutoff depending on the parameter specified
  if (parameter == "sample.size") {

    # Remove patients who dropped out of the trial or who are censored (they are not counted as completers)
    non.missing.outcome = !(is.na(selected.analysis.sample[, "outcome"]))
    completers = selected.analysis.sample[non.missing.outcome, ]

    # Check if there are any completers
    if (dim(completers)[1] > 0) {

      # Sort by the patient end time
      completers = completers[order(completers[, "patient.end.time"]), ]

      # Total number of completers
      total.sample.size = dim(completers)[1]

      # Find the time cutoff corresponding to the specified sample size in the selected sample
      # Truncate the VALUE argument if greater than the total number of completers
      time.cutoff = completers[min(value, total.sample.size), "patient.end.time"]

    } else {

      # No completers
      time.cutoff = 0
    }

  }

  if (parameter == "event") {

    # Remove patients who dropped out of the trial or who are censored (they are not counted as events)
    non.missing.outcome = !(is.na(selected.analysis.sample[, "outcome"])) & (selected.analysis.sample[, "patient.censor.indicator"] == 0)
    events = selected.analysis.sample[non.missing.outcome, ]

    # Check if there are any events
    if (dim(events)[1] > 0) {

      # Sort by the patient end time (this is when the events occurred)
      events = events[order(events[, "patient.end.time"]), ]

      # Total number of events
      total.event.count = dim(events)[1]

      # Find the time cutoff corresponding to the specified event count in the selected sample
      # Truncate the VALUE argument if greater than the total event count
      time.cutoff = events[min(value, total.event.count), "patient.end.time"]

    } else {

      # No events
      time.cutoff = 0
    }

  }

  if (parameter == "time") {

    # Time cutoff is directly specified
    time.cutoff = value

  }

  # Create the data slice by applying the time cutoff to all analysis samples

  # Loop over the analysis samples
  for (i in 1:n.analysis.samples) {

    sliced.analysis.sample = data.scenario$sample[[i]]$data

    if (data.scenario$sample[[i]]$outcome.type == "event") {

      # Apply slicing rules for event-type outcomes

      # If the patient end time is greater than the time cutoff, the patient censor indicator is set to TRUE
      sliced.analysis.sample[, "patient.censor.indicator"] = (sliced.analysis.sample[, "patient.end.time"] > time.cutoff)

      # If the patient end time or dropout time is greater than the time cutoff, the patient end time is set to the time cutoff
      sliced.analysis.sample[, "patient.end.time"] = pmin(time.cutoff, sliced.analysis.sample[, "patient.end.time"])
      sliced.analysis.sample[, "patient.dropout.time"] = pmin(time.cutoff, sliced.analysis.sample[, "patient.dropout.time"])

      # Outcome is truncated at the time cutoff for censored observations
      sliced.analysis.sample[, "outcome"] = (time.cutoff - sliced.analysis.sample[, "patient.start.time"]) * sliced.analysis.sample[, "patient.censor.indicator"] +
        sliced.analysis.sample[, "outcome"] * (1 - sliced.analysis.sample[, "patient.censor.indicator"])

      # If the patient outcome is negative, the outcome is set to NA (patient is enrolled after the time cutoff)
      sliced.analysis.sample[sliced.analysis.sample[, "outcome"] < 0, "outcome"] = NA

    } else {

      # Apply slicing rules for standard outcomes (binary and continuous outcome variables)

      # If the patient end time is greater than the time cutoff, the patient is considered a dropout and the outcome is set to NA
      sliced.analysis.sample[sliced.analysis.sample[, "patient.end.time"] > time.cutoff, "outcome"] = NA

      # If the patient end time or dropout time is greater than the time cutoff, the patient end time is set to the time cutoff
      sliced.analysis.sample[, "patient.end.time"] = pmin(time.cutoff, sliced.analysis.sample[, "patient.end.time"])
      sliced.analysis.sample[, "patient.dropout.time"] = pmin(time.cutoff, sliced.analysis.sample[, "patient.dropout.time"])

    }

    # Put the sliced data sample in the data slice
    data.slice$sample[[i]]$data = sliced.analysis.sample

  } # Loop over the analysis samples

  return(data.slice)
}
# End of CreateDataSlice
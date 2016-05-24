#######################################################################################################################

# Function: GeneratePatients.
# Argument: Design parameter, outcome parameter, sample id and number of patients or events to generate.
# Description: Generates data frames of simulated patients. This function is used in the CreateDataStack function.

GeneratePatients = function(current.design.parameter, current.outcome, current.sample.id, number){

  # Generate a set of outcome variables
  current.outcome.call = list(number, current.outcome$par)
  current.outcome.variables = as.matrix(do.call(current.outcome$dist, list(current.outcome.call)))
  colnames(current.outcome.variables) = paste0("outcome",1:ncol(current.outcome.variables))

  # Generate a set of design variables
  if (!is.null(current.design.parameter)){

    # Compute patient start times
    # Uniform patient start times
    if (current.design.parameter$enroll.dist == "UniformDist") {
      # Uniform distribution over [0, 1]
      enroll.par = list(number, list(max = 1))
      # Uniform distribution is expanded over the enrollment period
      patient.start.time = current.design.parameter$enroll.period * sort(unlist(lapply(list(enroll.par), "UniformDist")))
    } else {
      # Non-uniform patient start times
      # List of enrollment parameters
      enroll.par = list(number, current.design.parameter$enroll.dist.par)
      patient.start.time = sort(unlist(lapply(list(enroll.par), current.design.parameter$enroll.dist)))
    }
    # Patient start times are truncated at the end of the enrollment period
    patient.start.time = pmin(patient.start.time, current.design.parameter$enroll.period)

    # Compute patient end times
    # Patient end times
    if (!is.na(current.design.parameter$followup.period)) {
      # In a design with a fixed follow-up (followup.period is specified), the patient end time
      # is equal to the patient start time plus the fixed follow-up time
      patient.end.time = patient.start.time + current.design.parameter$followup.period
    }
    if (!is.na(current.design.parameter$study.duration)) {
      # In a design with a variable follow-up (study.duration is specified), the patient end time
      # is equal to the end of the trial
      patient.end.time = rep(current.design.parameter$study.duration, number)
    }

    # Compute patient dropout times (if the dropout distribution is specified) for the maximum sample size
    if (!is.na(current.design.parameter$dropout.dist)) {
      # Uniform patient dropout times
      if (current.design.parameter$dropout.dist == "UniformDist") {
        # Uniform distribution over [0, 1]
        dropout.par = list(number, 1/current.design.parameter$dropout.dist.par)
        # Uniform distribution is expanded over the patient-specific periods
        patient.dropout.time = patient.start.time + (patient.end.time - patient.start.time) *
          unlist(lapply(list(dropout.par), "UniformDist"))
      } else {
        # Non-uniform patient dropout times
        # List of dropout parameters
        dropout.par = list(number, current.design.parameter$dropout.dist.par)
        patient.dropout.time = patient.start.time +
          unlist(lapply(list(dropout.par), current.design.parameter$dropout.dist))
      }
      # If the patient end time is greater than the patient dropout time, the patient end time
      # is truncated, the patient dropout indicator is set to TRUE.
      patient.dropout.indicator = (patient.end.time >= patient.dropout.time)
      patient.end.time = pmin(patient.end.time, patient.dropout.time)
    } else {
      # No dropout distribution is specified
      patient.dropout.time = rep(NA, number)
      patient.dropout.indicator = rep(FALSE, number)
    }

    # Patient censore will be get later on in the function according to the outcome variable
    patient.censor.indicator = rep(FALSE, number)

    # Create a data frame and save it
    current.design.variables = t(rbind(patient.start.time,
                                       patient.end.time,
                                       patient.dropout.time,
                                       patient.dropout.indicator,
                                       patient.censor.indicator))

  } else if (is.null(current.design.parameter)){
    # No design parameters are specified in the data model
    patient.start.time = rep(NA, number)
    patient.end.time = rep(NA, number)
    patient.dropout.time = rep(NA, number)
    patient.dropout.indicator = rep(FALSE, number)
    patient.censor.indicator = rep(FALSE, number)

    # Create a data frame and save it
    current.design.variables = t(rbind(patient.start.time,
                                       patient.end.time,
                                       patient.dropout.time,
                                       patient.dropout.indicator,
                                       patient.censor.indicator))

  }

  colnames(current.design.variables) = c("patient.start.time",
                                         "patient.end.time",
                                         "patient.dropout.time",
                                         "patient.dropout.indicator",
                                         "patient.censor.indicator")

  # Create the list with the data frame for the current design and outcome parameter and for each outcome
  current.design.outcome.variables = list()

  # Create the censor indicator for each outcome
  for (outcome.index in 1:length(current.outcome$type)){

    current.outcome.type = current.outcome$type[outcome.index]
    patient.end.time = current.design.variables[,"patient.end.time"]
    patient.start.time = current.design.variables[,"patient.start.time"]
    patient.dropout.time = current.design.variables[,"patient.dropout.time"]
    patient.censor.indicator = current.design.variables[,"patient.censor.indicator"]
    outcome = current.outcome.variables[,paste0("outcome",outcome.index)]

    # Compute patient censor times for the analysis data sample if the current outcome type is "event"
    if (current.outcome.type == "event") {

      # Dropout distribution is specified
      if (!all(is.na(patient.dropout.time))) {

        # Outcome variable is truncated and the patient censor indicator is set to TRUE
        # if the outcome variable is greater than the patient dropout time (relative to the patient start time)
        outcome = pmin(outcome, patient.dropout.time - patient.start.time)

      }

      # Enrollment distribution is specified
      if (!all(is.na(patient.start.time))) {

        patient.censor.indicator = patient.dropout.indicator

        # Outcome variable is truncated and the patient censor indicator is set to TRUE
        # if the outcome variable is greater than the patient end time (relative to the patient start time)
        patient.censor.indicator = patient.censor.indicator | (outcome >= patient.end.time - patient.start.time)
        outcome = pmin(outcome, patient.end.time - patient.start.time)

        # Patient end time (relative to the patient start time) is set to the outcome variable if the
        # patient experience the event (that is, the patient censor indicator is FALSE)
        patient.end.time = (!patient.censor.indicator) * (patient.start.time + outcome) +
          (patient.censor.indicator) * patient.end.time
      }

    } else {
      # Current outcome type is "standard"

      # Dropout distribution is specified
      if (!all(is.na(patient.dropout.time))) {

        # Outcome variable is set to NA if the patient dropout indicator is TRUE
        outcome[patient.dropout.indicator] = NA
      }

      patient.censor.indicator = rep(FALSE, length(patient.censor.indicator))
    }

    # Create a data frame for the current sample and outcome
    df = t(rbind(outcome,
                 patient.start.time,
                 patient.end.time,
                 patient.dropout.time,
                 patient.censor.indicator))

    colnames(df) = c("outcome",
                     "patient.start.time",
                     "patient.end.time",
                     "patient.dropout.time",
                     "patient.censor.indicator")

    current.design.outcome.variables[[outcome.index]] = list(id = current.sample.id[outcome.index],
                                                             outcome.type = current.outcome.type,
                                                             data = df)

  }

  return(current.design.outcome.variables)
} # End of GeneratePatients function

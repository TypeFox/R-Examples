######################################################################################################################

# Function: CreateDataStructure.
# Argument: Data model.
# Description: This function is based on the old data_model_extract function. It performs error checks in the data model
# and creates a "data structure", which is an internal representation of the original data model used by all other Mediana functions.

CreateDataStructure = function(data.model) {

  # Check the general set
  if (is.null(data.model$samples))
    stop("Data model: At least one sample must be specified.")

  # Number of samples in the data model
  n.samples = length(data.model$samples)

  if (is.null(data.model$general))
    stop("Data model: General set of parameters must be specified.")

  # General set of parameters

  # List of outcome distribution parameters
  outcome = list()

  # Outcome distribution is required in the general set of data model parameters
  if (is.null(data.model$general$outcome.dist))
    stop("Data model: Outcome distribution must be specified in the general set of parameters.")
  outcome.dist = data.model$general$outcome.dist
  if (!exists(outcome.dist)) {
    stop(paste0("Data model: Outcome distribution function '", outcome.dist, "' does not exist."))
  } else {
    if (!is.function(get(as.character(outcome.dist), mode = "any")))
      stop(paste0("Data model: Outcome distribution function '", outcome.dist, "' does not exist."))
  }

  # Extract sample-specific parameters

  # List of outcome parameter sets
  outcome.parameter.set = list()
  # List of design parameter sets
  design.parameter.set = list()
  # List of sample IDs
  id = list()

  # Determine if the data model is expanded or compact (compact if the sample size sets are
  # specified in the general set of parameters, extended if the sample size sets
  # are specified for each sample)
  compact.size = FALSE
  expanded.size = FALSE

  sample.size = FALSE
  event = FALSE

  if (is.null(data.model$general$sample.size) & is.null(data.model$general$event)) {
    if (is.null(data.model$samples[[1]]$sample.size) & is.null(data.model$samples[[1]]$event))
      stop("Data model: Sample sizes or events must be specified either in the general set or in the sample-specific set of parameters.")
  }
  if (!is.null(data.model$general$sample.size)) {
    if (!is.null(data.model$samples[[1]]$sample.size))
      stop("Data model: Sample sizes must be specified either in the general set or in the sample-specific set of parameters but not both.")
  }
  if (!is.null(data.model$general$event)) {
    if (!is.null(data.model$samples[[1]]$event))
      stop("Data model: Events must be specified either in the general set or in the sample-specific set of parameters but not both.")
  }

  if (!is.null(data.model$general$event) & !is.null(data.model$general$sample.size)) {
    stop("Data model: Sample sizes or Events must be specified but not both.")
  }

  if (!is.null(data.model$samples[[1]]$event) & !is.null(data.model$samples[[1]]$sample.size)) {
    stop("Data model: Sample sizes or Events must be specified but not both.")
  }

  if (!is.null(data.model$samples[[1]]$event) & !is.null(data.model$general$sample.size)) {
    stop("Data model: Sample sizes or Events must be specified but not both.")
  }

  if (!is.null(data.model$general$event) & !is.null(data.model$samples[[1]]$sample.size)) {
    stop("Data model: Sample sizes or Events must be specified but not both.")
  }


  # Compute the number of sample size sets
  if (!is.null(data.model$general$sample.size) | !is.null(data.model$samples[[1]]$sample.size)){
    sample.size = TRUE
    if (!is.null(data.model$general$sample.size)) {
      compact.size = TRUE
      n.sample.size.sets = length(data.model$general$sample.size)
    } else {
      expanded.size = TRUE
      n.sample.size.sets = length(data.model$samples[[1]]$sample.size)
      for (i in 1:n.samples) {
        if (is.null(data.model$samples[[i]]$sample.size))
          stop("Data model: Sample sizes must be specified for all samples.")
        if (n.sample.size.sets != length(data.model$samples[[i]]$sample.size))
          stop("Data model: The same number of sample sizes must be specified across the samples.")
      }
    }

    # Data frame of sample size sets
    sample.size.set = matrix(0, n.sample.size.sets, n.samples)

    # Create a list of sample size sets
    for (i in 1:n.sample.size.sets) {
      if (expanded.size) {
        for (j in 1:n.samples) {
          sample.size.set[i, j] = data.model$samples[[j]]$sample.size[[i]]
        }
      }
      if (compact.size) {
        for (j in 1:n.samples) {
          sample.size.set[i, j] = data.model$general$sample.size[[i]]
        }
      }
    }
    sample.size.set = as.data.frame(sample.size.set)

    # Error check
    if (any(sample.size.set<=0)) stop("Data model : Sample size must be strictly positive")

  } else  {
    sample.size.set = NA
  }

  # Compute the number of event sets
  if (!is.null(data.model$general$event)){
    event = TRUE
    compact.size = TRUE
    event.set = data.frame(event.total = data.model$general$event$n.events)
    rando.ratio = data.model$general$event$rando.ratio
    if (is.null(rando.ratio)) rando.ratio = rep(1,n.samples)

    # Error check
    if (any(event.set<=0)) stop("Data model : Number of events must be strictly positive")
    if (length(rando.ratio) != n.samples) stop("Data model: the randomization ratio of each sample must be specified")
    if (any(rando.ratio<=0)) stop("Data model: the randomization ratio of each sample must be positive")
    if (any(rando.ratio %%1 != 0)) stop("Data model: the randomization ratio of each sample must be an integer")

  } else  {
    event.set = NA
    rando.ratio = NA
  }

  # Compute the number of outcome parameter sets
  for (i in 1:n.samples) {
    if (is.null(data.model$samples[[i]]$outcome.par))
      stop("Data model: Outcome parameters must be specified for all samples.")

    outcome.par = data.model$samples[[i]]$outcome.par

    if (i == 1) {
      n.outcome.parameter.sets = length(outcome.par)
    } else {
      if (n.outcome.parameter.sets != length(outcome.par))
        stop("Data model: The same number of outcome parameter sets must be specified across the samples.")
    }
  }

  # Create a list of outcome parameter sets
  for (i in 1:n.outcome.parameter.sets) {
    temp = list()
    for (j in 1:n.samples) {
      temp[[j]] = data.model$samples[[j]]$outcome.par[[i]]
      # Check if the outcome parameters are correctly specified and determine the dimensionality of the outcome distribution
      dummy.function.call = list(1, data.model$samples[[j]]$outcome.par[[i]])
      outcome.dist.dim = length(do.call(outcome.dist, list(dummy.function.call)))
    }
    outcome.parameter.set[[i]] = temp
  }

  if (is.null(data.model$general$outcome.type) & sample.size == TRUE) {
    outcome.type = rep("standard", outcome.dist.dim)
  } else if (is.null(data.model$general$outcome.type) & event == TRUE) {
    outcome.type = rep("event", outcome.dist.dim)
  } else {
    outcome.type = data.model$general$outcome.type
    if (length(outcome.type) != outcome.dist.dim)
      stop("Data model: Number of outcome types must be equal to the number of dimensions in the outcome distribution.")
  }

  # Create a list of sample IDs
  for (i in 1:n.samples) {
    if (is.null(data.model$samples[[i]]$id))
      stop("Data model: Sample IDs must be specified for all samples.")
    if (outcome.dist.dim != length(data.model$samples[[i]]$id))
      stop("Data model: The same number of sample IDs in each sample must be equal to the number of dimensions in the outcome distribution.")
    id[[i]] = data.model$samples[[i]]$id
  }

  # Compute the number of design parameter sets
  if (is.null(data.model$general$design)) {
    n.design.parameter.sets = NA
    design.parameter.set = NULL
  } else {
    n.design.parameter.sets = length(data.model$general$design)
  }

  # Create a list of design parameter sets
  if (!is.null(design.parameter.set)) {
    for (i in 1:n.design.parameter.sets) {
      if (!is.null(data.model$general$design[[i]]$followup.period) & !is.null(data.model$general$design[[i]]$study.duration))
        stop("Data model: Either the length of the follow-up period or total study duration can be specified but not both.")

      if (is.null(data.model$general$design[[i]]$enroll.dist) & !is.null(data.model$general$design[[i]]$dropout.dist))
        stop("Data model: Dropout parameters may not be specified without enrollment parameters.")

      if (is.null(data.model$general$design[[i]]$enroll.period)) {
        enroll.period = NA
      } else {
        enroll.period = data.model$general$design[[i]]$enroll.period
      }

      if (is.null(data.model$general$design[[i]]$enroll.dist)) {
        enroll.dist = NA
      } else {
        enroll.dist = data.model$general$design[[i]]$enroll.dist
        if (!exists(enroll.dist)) {
          stop(paste0("Data model: Enrollment distribution function '", enroll.dist, "' does not exist."))
        } else {
          if (!is.function(get(as.character(enroll.dist), mode = "any")))
            stop(paste0("Data model: Enrollment distribution function '", enroll.dist, "' does not exist."))
        }
      }

      if (enroll.dist == "UniformDist") {
        enroll.dist.par = NA
      } else {
        if (is.null(data.model$general$design[[i]]$enroll.dist.par)) {
          stop("Data model: Enrollment distribution parameters must be specified for non-uniform distributions.")
        } else {
          enroll.dist.par = data.model$general$design[[i]]$enroll.dist.par
        }
      }

      if (is.null(data.model$general$design[[i]]$followup.period)) {
        followup.period = NA
      } else {
        followup.period = data.model$general$design[[i]]$followup.period
      }

      if (is.null(data.model$general$design[[i]]$study.duration)) {
        study.duration = NA
      } else {
        study.duration = data.model$general$design[[i]]$study.duration
      }

      if (is.null(data.model$general$design[[i]]$dropout.dist)) {
        dropout.dist = NA
      } else {
        dropout.dist = data.model$general$design[[i]]$dropout.dist
        if (!exists(dropout.dist)) {
          stop(paste0("Data model: Dropout distribution function '", dropout.dist, "' does not exist."))
        } else {
          if (!is.function(get(as.character(dropout.dist), mode = "any")))
            stop(paste0("Data model: Dropout distribution function '", dropout.dist, "' does not exist."))
        }
      }

      if (is.null(data.model$general$design[[i]]$dropout.dist.par)) {
        dropout.dist.par = NA
      } else {
        dropout.dist.par = data.model$general$design[[i]]$dropout.dist.par
      }

      design.parameter.set[[i]] = list(enroll.period = enroll.period,
                                       enroll.dist = enroll.dist,
                                       enroll.dist.par = enroll.dist.par,
                                       followup.period = followup.period,
                                       study.duration = study.duration,
                                       dropout.dist = dropout.dist,
                                       dropout.dist.par = dropout.dist.par)
    }
  }

  # Create the data structure
  outcome = list(outcome.dist = outcome.dist, outcome.type = outcome.type, outcome.dist.dim = outcome.dist.dim)
  data.structure = list(description = "data.structure",
                        id = id,
                        outcome = outcome,
                        sample.size.set = sample.size.set,
                        event.set = event.set,
                        rando.ratio = rando.ratio,
                        outcome.parameter.set = outcome.parameter.set,
                        design.parameter.set = design.parameter.set)
  return(data.structure)
}
# End of CreateDataStructure

######################################################################################################################

# Function: CreateAnalysisStructure.
# Argument: Analysis model.
# Description: This function is based on the old analysis_model_extract function. It performs error checks in the analysis model
# and creates an "analysis structure", which is an internal representation of the original analysis model used by all other Mediana functions.

CreateAnalysisStructure = function(analysis.model) {

  # Check the general set
  if (is.null(analysis.model$tests) & is.null(analysis.model$statistics))
    stop("Analysis model: At least one test or statistic must be specified.")

  # General set of analysis model parameters

  # Extract interim analysis parameters
  if (!is.null(analysis.model$general$interim.analysis)) {

    interim.looks = analysis.model$general$interim.analysis$interim.looks
    if (!(interim.looks$parameter %in% c("sample.size", "event", "time")))
      stop("Analysis model: Parameter in the interim analysis specifications must be sample.size, event or time.")
    interim.analysis = list(interim.looks = interim.looks)

  } else {
    interim.analysis = NULL
  }

  # Extract test-specific parameters

  if (!is.null(analysis.model$tests)) {

    # Number of tests in the analysis model
    n.tests = length(analysis.model$tests)

    # List of tests (id, statistical method, sample list, parameters)
    test = list()

    for (i in 1:n.tests) {
      # Test IDs
      if (is.null(analysis.model$tests[[i]]$id))
        stop("Analysis model: IDs must be specified for all tests.") else id = analysis.model$tests[[i]]$id
        # List of samples
        if (is.null(analysis.model$tests[[i]]$samples))
          stop("Analysis model: Samples must be specified for all tests.") else samples = analysis.model$tests[[i]]$samples
          # Statistical method
          method = analysis.model$test[[i]]$method

          if (!exists(method)) {
            stop(paste0("Analysis model: Statistical method function '", method, "' does not exist."))
          } else if (!is.function(get(as.character(method), mode = "any"))) {
            stop(paste0("Analysis model: Statistical method function '", method, "' does not exist."))
          }

          # Test parameters (optional)
          if (is.null(analysis.model$tests[[i]]$par)) par = NA else par = analysis.model$tests[[i]]$par

          test[[i]] = list(id = id, method = method, samples = samples, par = par)
    }

    # Check if id is uniquely defined
    if (any(table(unlist(lapply(test,function(list) list$id)))>1))
      stop("Analysis model: Tests IDs must be uniquely defined.")

  } else {
    # No tests are specified
    test = NULL

  }

  # Extract statistic-specific parameters

  if (!is.null(analysis.model$statistics)) {

    # Number of statistics in the analysis model
    n.statistics = length(analysis.model$statistics)

    # List of statistics (id, statistical method, sample list, parameters)
    statistic = list()

    for (i in 1:n.statistics) {
      # Statistic IDs
      if (is.null(analysis.model$statistic[[i]]$id))
        stop("Analysis model: IDs must be specified for all statistics.") else id = analysis.model$statistic[[i]]$id
        # List of samples
        if (is.null(analysis.model$statistic[[i]]$samples))
          stop("Analysis model: Samples must be specified for all statistics.") else samples = analysis.model$statistic[[i]]$samples
          # Statistical method
          method = analysis.model$statistic[[i]]$method

          if (!exists(method)) {
            stop(paste0("Analysis model: Statistical method function '", method, "' does not exist."))
          } else if (!is.function(get(as.character(method), mode = "any"))) {
            stop(paste0("Analysis model: Statistical method function '", method, "' does not exist."))
          }

          if (is.null(analysis.model$statistic[[i]]$par)) par = NA else par = analysis.model$statistic[[i]]$par

          statistic[[i]] = list(id = id, method = method, samples = samples, par = par)
    }

  } else {
    # No statistics are specified
    statistic = NULL

  }

  # Extract parameters of multiplicity adjustment methods

  # List of multiplicity adjustments (procedure, parameters, tests)
  mult.adjust = list(list())

  # Number of multiplicity adjustment methods
  if (is.null(analysis.model$general$mult.adjust)) {
    # No multiplicity adjustment is specified
    mult.adjust = NULL
  } else {
    n.mult.adjust = length(analysis.model$general$mult.adjust)
    for (i in 1:n.mult.adjust) {
      mult.adjust.temp = list()
      # Number of multiplicity adjustments within each mult.adj scenario
      n.mult.adj.sc=length(analysis.model$general$mult.adjust[[i]])
      for (j in 1:n.mult.adj.sc){
        proc = analysis.model$general$mult.adjust[[i]][[j]]$proc
        if (is.na(proc) | is.null(analysis.model$general$mult.adjust[[i]][[j]]$par)) par = NA else par = analysis.model$general$mult.adjust[[i]][[j]]$par
        if (is.null(analysis.model$general$mult.adjust[[i]][[j]]$tests)) {
          tests = lapply(test, function(list) list$id)
        } else {
          tests = analysis.model$general$mult.adjust[[i]][[j]]$tests
        }

        # If the multiplicity adjustment procedure is specified, check if it exists
        if (!is.na(proc)) {
          if (!exists(proc)) {
            stop(paste0("Analysis model: Multiplicity adjustment procedure function '", proc, "' does not exist."))
          } else if (!is.function(get(as.character(proc), mode = "any"))) {
            stop(paste0("Analysis model: Multiplicity adjustment procedure function '", proc, "' does not exist."))
          }
        }

        # Check if tests defined in the multiplicity adjustment exist (defined in the test list)
        temp_list = lapply(lapply(tests,function(l1,l2) l1 %in% l2, lapply(test, function(list) list$id)), function(l) any(l == FALSE))

        if (!is.na(proc) & any(temp_list == TRUE))
          stop(paste0("Analysis model: Multiplicity adjustment procedure test has not been specified in the test-specific model."))

        mult.adjust.temp[[j]] = list(proc = proc, par = par, tests = tests)
      }

      mult.adjust[[i]] = mult.adjust.temp
      # Check if tests defined in multiplicity adjustment is defined in one and only one multiplicity adjustment
      if (any(table(unlist(lapply(mult.adjust[[i]],function(list) list$tests)))>1))
        stop(paste0("Analysis model: Multiplicity adjustment procedure test has been specified in more than one multiplicity adjustment."))

    }
  }

  # Create the analysis structure
  analysis.structure = list(description = "analysis.structure",
                            test = test,
                            statistic = statistic,
                            mult.adjust = mult.adjust,
                            interim.analysis = interim.analysis)
  return(analysis.structure)

}
# End of CreateAnalysisStructure
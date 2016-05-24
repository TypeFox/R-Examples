############################################################################################################################

# Function: CSE.default
# Argument: Data model (or data stack), analysis model (or analysis stack) and evaluation model.
# Description: This function applies the metrics specified in the evaluation model to the test results (p-values) and
# summaries to the statistic results.
#' @export
CSE.default = function(data, analysis, evaluation, simulation) {

  # Error check
  if (!(class(data) %in% c("DataModel", "DataStack"))) stop("CSE: a DataModel object must be specified in the data parameter")
  if (!(class(analysis) %in% c("AnalysisModel", "AnalysisStack"))) stop("CSE: an AnalysisModel object must be specified in the analysis parameter")
  if (!(class(evaluation) %in% c("EvaluationModel"))) stop("CSE: an EvaluationModel object must be specified in the evaluation parameter")
  if (!(class(simulation) %in% c("SimParameters"))) stop("CSE: a SimParameters object must be specified in the simulation parameter")

  # Start time
  start.time = Sys.time()


  # Perform error checks for the evaluation model and create an internal evaluation structure
  #	(in the first position in order to be sure that if any error is made, the simulation won't run)
  evaluation.structure = CreateEvaluationStructure(evaluation)

  # Case 1: Data model and Analysis model
  if (class(data) == "DataModel" & class(analysis) == "AnalysisModel"){
    data.model = data
    analysis.model = analysis

    # Data structure
    data.structure = CreateDataStructure(data.model)

    # Create the analysis stack from the specified data and analysis models
    analysis.stack = PerformAnalysis(data.model, analysis.model, sim.parameters = simulation)

    # Analysis structure
    analysis.structure = analysis.stack$analysis.structure

    # Simulation parameters
    sim.parameters = analysis.stack$sim.parameters
  }

  # Case 2: Data stack and Analysis model
  if (class(data) == "DataStack" & class(analysis) == "AnalysisModel"){
    data.stack = data
    analysis.model = analysis

    # Data structure
    data.structure = data.stack$data.structure

    # Create the analysis stack from the specified data and analysis models
    analysis.stack = PerformAnalysis(data.stack, analysis.model, sim.parameters = simulation)

    # Analysis structure
    analysis.structure = analysis.stack$analysis.structure

    # Simulation parameters
    sim.parameters = analysis.stack$sim.parameters
  }

  # Case 3: Data stack and Analysis stack
  if (class(data) == "DataStack" & class(analysis) == "AnalysisStack"){
    data.stack = data
    analysis.stack = analysis

    # Data structure
    data.structure = data.stack$data.structure

    # Analysis structure
    analysis.structure = analysis.stack$analysis.structure

    # Simulation parameters
    if (!is.null(simulation))
      warning("The simulation parameters (simulation) defined in the CSE function will be ignored as an analysis stack has been defined.")
    sim.parameters = analysis.stack$sim.parameters
  }

  # Get the number of simulation runs
  n.sims = sim.parameters$n.sims


  # Extract the analysis scenario grid and compute the number of analysis scenarios in this analysis stack
  analysis.scenario.grid = analysis.stack$analysis.scenario.grid
  n.analysis.scenarios = dim(analysis.scenario.grid)[1]

  # Number of outcome parameter sets, sample size sets and design parameter sets
  n.outcome.parameter.sets = length(levels(as.factor(analysis.scenario.grid$outcome.parameter)))
  n.design.parameter.sets = length(levels(as.factor(analysis.scenario.grid$design.parameter)))
  n.sample.size.sets = length(levels(as.factor(analysis.scenario.grid$sample.size)))

  # Number of data scenario
  n.data.scenarios =  n.outcome.parameter.sets*n.sample.size.sets*n.design.parameter.sets

  # Number of multiplicity adjustment procedure
  n.mult.adj = length(levels(as.factor(analysis.scenario.grid$multiplicity.adjustment)))

  # Number of criteria specified in the evaluation model
  n.criteria = length(evaluation.structure$criterion)

  # Criterion IDs
  criterion.id = rep(" ", n.criteria)

  # Check if the tests and statistics referenced in the evaluation model are actually defined in the analysis model

  # Number of tests specified in the analysis model
  n.tests = length(analysis.structure$test)

  # Number of statistics specified in the analysis model
  n.statistics = length(analysis.structure$statistic)

  if(is.null(analysis.structure$test)) {

    # Test IDs
    test.id = " "

  } else {

    # Test IDs
    test.id = rep(" ", n.tests)

    for (test.index in 1:n.tests) {
      test.id[test.index] = analysis.structure$test[[test.index]]$id
    }

  }

  if(is.null(analysis.structure$statistic)) {

    # Statistic IDs
    statistic.id = " "

  } else {

    # Statistic IDs
    statistic.id = rep(" ", n.statistics)

    for (statistic.index in 1:n.statistics) {
      statistic.id[statistic.index] = analysis.structure$statistic[[statistic.index]]$id
    }

  }

  for (criterion.index in 1:n.criteria) {

    current.criterion = evaluation.structure$criterion[[criterion.index]]

    criterion.id[criterion.index] = current.criterion$id

    # Number of tests specified within the current criterion
    n.criterion.tests = length(current.criterion$tests)

    # Number of statistics specified within the current criterion
    n.criterion.statistics = length(current.criterion$statistics)

    if (n.criterion.tests > 0) {
      for (i in 1:n.criterion.tests) {
        if (!(current.criterion$tests[[i]] %in% test.id))
          stop(paste0("Evaluation model: Test '", current.criterion$tests[[i]], "' is not defined in the analysis model."))
      }
    }

    if (n.criterion.statistics > 0) {
      for (i in 1:n.criterion.statistics) {
        if (!(current.criterion$statistics[[i]] %in% statistic.id))
          stop(paste0("Evaluation model: Statistic '", current.criterion$statistics[[i]], "' is not defined in the analysis model."))
      }
    }

  }

  # Number of analysis points (total number of interim and final analyses)
  if (!is.null(analysis.structure$interim.analysis)) {
    n.analysis.points = length(analysis.structure$interim.analysis$interim.looks$fraction)
  } else {
    # No interim analyses
    n.analysis.points = 1
  }

  # Create the evaluation stack (list of evaluation results for each analysis scenario in the analysis stack)
  #evaluation.set = list()

  # List of analysis scenarios
  analysis.scenario = list()

  analysis.scenario.index = 0

  # Loop over the data scenarios
  for (data.scenario.index in 1:n.data.scenarios) {

    # Loop over the multiplicity adjustment
    for (mult.adj.index in 1:n.mult.adj) {

      analysis.scenario.index =  analysis.scenario.index +1

      # List of criteria
      criterion = list()


      # Loop over the criteria
      for (criterion.index in 1:n.criteria) {


        # Current metric
        current.criterion = evaluation.structure$criterion[[criterion.index]]

        # Number of tests specified in the current criterion
        n.criterion.tests = length(current.criterion$tests)

        # Number of statistics specified in the current criterion
        n.criterion.statistics = length(current.criterion$statistics)

        # Create a matrix of test results (p-values) across the simulation runs for the current criterion and analysis scenario
        if (n.criterion.tests > 0) {

          test.result.matrix = matrix(0, n.sims, n.criterion.tests)

          # Create a template for selecting the test results (p-values)
          test.result.flag = lapply(analysis.structure$test, function(x) any(current.criterion$tests == x$id))
          test.result.flag.num = match(current.criterion$tests,test.id)


        } else {
          test.result.matrix = NA
          test.result.flag = NA
          test.result.flag.num = NA
        }

        # Create a matrix of statistic results across the simulation runs for the current criterion and analysis scenario
        if (n.criterion.statistics > 0) {

          statistic.result.matrix = matrix(0, n.sims, n.criterion.statistics)

          # Create a template for selecting the statistic results
          statistic.result.flag = lapply(analysis.structure$statistic, function(x) any(current.criterion$statistics == x$id))
          statistic.result.flag.num = match(current.criterion$statistics ,statistic.id)

        } else {
          statistic.result.matrix = NA
          statistic.result.flag = NA
          statistic.result.flag.num = NA
        }

        # Loop over the analysis points
        for(analysis.point.index in 1:n.analysis.points){

          # Loop over the simulation runs
          for (sim.index in 1:n.sims) {

            # Current analysis scenario
            current.analysis.scenario = analysis.stack$analysis.set[[sim.index]][[data.scenario.index]]$result$tests.adjust$analysis.scenario[[mult.adj.index]]

            # Extract the tests specified in the current criterion
            if (n.criterion.tests > 0) {
              #test.result.matrix[sim.index, ] = current.analysis.scenario[test.result.flag, analysis.point.index]
              test.result.matrix[sim.index, ] = current.analysis.scenario[test.result.flag.num, analysis.point.index]
            }

            # Extract the statistics specified in the current criterion
            if (n.criterion.statistics > 0) {
              #statistic.result.matrix[sim.index, ] = analysis.stack$analysis.set[[sim.index]][[data.scenario.index]]$result$statistic[statistic.result.flag, analysis.point.index]
              statistic.result.matrix[sim.index, ] = analysis.stack$analysis.set[[sim.index]][[data.scenario.index]]$result$statistic[statistic.result.flag.num, analysis.point.index]
            }

          } # Loop over the simulation runs

          # Apply the method specified in the current metric with metric parameters
          single.result = as.matrix(do.call(current.criterion$method,
                                            list(test.result.matrix, statistic.result.matrix, current.criterion$par)))

          if (n.analysis.points == 1) {
            # Only one analysis point (final analysis) is specified
            evaluation.results = single.result
          } else {
            # Two or more analysis points (interim and final analyses) are specified
            if (analysis.point.index == 1) {
              evaluation.results = single.result
            } else {
              evaluation.results = cbind(evaluation.results, single.result)
            }
          }

        } # Loop over the analysis points

        evaluation.results = as.data.frame(evaluation.results)

        # Assign labels
        rownames(evaluation.results) = unlist(current.criterion$labels)

        if (n.analysis.points == 1) {
          colnames(evaluation.results) = "Analysis"
        } else {
          names = rep("", n.analysis.points)
          for (j in 1:n.analysis.points) names[j] = paste0("Analysis ", j)
          colnames(evaluation.results) = names
        }

        criterion[[criterion.index]] = list(id = criterion.id[[criterion.index]],
                                            result = evaluation.results)

      } # Loop over the criteria

      analysis.scenario[[analysis.scenario.index]] = list(criterion = criterion)

    } # Loop over the multiplicity adjustment

  } # Loop over the data scenarios

  #evaluation.set = list(analysis.scenario = analysis.scenario)

  # Create a single data frame with simulation results
  simulation.results = data.frame(sample.size = numeric(),
                                  outcome.parameter = numeric(),
                                  design.parameter = numeric(),
                                  multiplicity.adjustment = numeric(),
                                  criterion = character(),
                                  test.statistic = character(),
                                  result = numeric(),
                                  stringsAsFactors = FALSE)

  row.index = 1

  n.analysis.scenarios = length(analysis.scenario)

  for (scenario.index in 1:(n.data.scenarios*n.mult.adj)) {

    current.analysis.scenario = analysis.scenario[[scenario.index]]
    n.criteria = length(current.analysis.scenario$criterion)
    current.analysis.scenario.grid = analysis.scenario.grid[scenario.index, ]

    for (j in 1:n.criteria) {
      n.rows = dim(current.analysis.scenario$criterion[[j]]$result)[1]
      for (k in 1:n.rows) {
        simulation.results[row.index, 1] = current.analysis.scenario.grid[1, 3]
        simulation.results[row.index, 2] = current.analysis.scenario.grid[1, 2]
        simulation.results[row.index, 3] = current.analysis.scenario.grid[1, 1]
        simulation.results[row.index, 4] = current.analysis.scenario.grid[1, 4]
        simulation.results[row.index, 5] = current.analysis.scenario$criterion[[j]]$id
        simulation.results[row.index, 6] = rownames(current.analysis.scenario$criterion[[j]]$result)[k]
        simulation.results[row.index, 7] = current.analysis.scenario$criterion[[j]]$result[k, 1]
        row.index = row.index + 1
      }
    }


  }

  end.time = Sys.time()
  timestamp  = list(start.time = start.time,
                    end.time = end.time,
                    duration = difftime(end.time,start.time, units = "mins"))

  # Create the evaluation stack
  evaluation.stack = list(#description = "CSE",
    simulation.results = simulation.results,
    #evaluation.set = evaluation.set,
    analysis.scenario.grid = analysis.scenario.grid,
    data.structure = data.structure,
    analysis.structure = analysis.structure,
    evaluation.structure = evaluation.structure,
    sim.parameters = sim.parameters,
    #env.information = env.information,
    timestamp = timestamp )

  class(evaluation.stack) = "CSE"
  return(evaluation.stack)

}
# End of CSE
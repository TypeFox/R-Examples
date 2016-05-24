## niching.R
##   - Functions defining evolution main loops for cluster-based niching GP
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Cluster-based multi-niche genetic programming
##'
##' Perform a multi-niche genetic programming run. The required argument
##' \code{fitnessFunction} must be supplied with an objective function that assigns
##' a numerical fitness value to an R function. Fitness values are minimized, i.e.
##' smaller values mean higher/better fitness. If a multi-objective
##' \code{selectionFunction} is used, \code{fitnessFunction} return a numerical
##' vector of fitness values.
##' In a multi-niche genetic programming run, the initial population is clustered
##' via a \code{clusterFunction} into \code{numberOfNiches} niches. In each niche,
##' a genetic programming run is executed with \code{passStopCondition} as stop
##' condition. These runs are referred to as a parallel pass. After each parallel
##' pass, the niches are joined again using a \code{joinFunction} into a population.
##' From here, the process starts again with a clustering step, until the global
##' \code{stopCondition} is met.
##' The result of the multi-niche genetic programming run is a genetic programming
##' result object containing a GP population of R functions.
##' \code{summary.geneticProgrammingResult} can be used to create summary views of a
##' GP result object.
##'
##' @param fitnessFunction In case of a single-objective selection function,
##'   \code{fitnessFunction} must be a single function that assigns a
##'   numerical fitness value to a GP individual represented as a R function.
##'   Smaller fitness values mean higher/better fitness. If a multi-objective
##'   selection function is used, \code{fitnessFunction} must return a numerical
##'   vector of fitness values.
##' @param stopCondition The stop condition for the evolution main loop. See
##'   \code{makeStepsStopCondition} for details.
##' @param passStopCondition The stop condition for each parallel pass. See
##'   \code{makeStepsStopCondition} for details.
##' @param numberOfNiches The number of niches to cluster the population into.
##' @param clusterFunction The function used to cluster the population into
##'   niches. The first parameter of this function is a GP population, the
##'   second paramater an integer representing the number of niches. Defaults
##'   to \code{\link{groupListConsecutive}}.
##' @param joinFunction The function used to join all niches into a population
##'   again after a round of parallel passes. Defaults to a function that
##'   simply concatenates all niches.
##' @param population The GP population to start the run with. If this parameter
##'   is missing, a new GP population of size \code{populationSize} is created
##'   through random growth.
##' @param populationSize The number of individuals if a population is to be
##'   created.
##' @param eliteSize The number of "elite" individuals to keep. Defaults to
##'  \code{ceiling(0.1 * populationSize)}.
##' @param elite The elite list, must be alist of individuals sorted in ascending
##'   order by their first fitness component.
##' @param functionSet The function set.
##' @param inputVariables The input variable set.
##' @param constantSet The set of constant factory functions.
##' @param searchHeuristic The search-heuristic (i.e. optimization algorithm) to use
##'   in the search of solutions. See the documentation for \code{searchHeuristics} for
##'   available algorithms.
##' @param crossoverFunction The crossover function.
##' @param mutationFunction The mutation function.
##' @param restartCondition The restart condition for the evolution main loop. See
##'   \link{makeFitnessStagnationRestartCondition} for details.
##' @param restartStrategy The strategy for doing restarts. See
##'   \link{makeLocalRestartStrategy} for details.
##' @param progressMonitor A function of signature
##'   \code{function(population, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber,
##'   bestFitness, timeElapsed, ...)} to be called with each evolution step. Seach heuristics
##'   may pass additional information via the \code{...} parameter.
##' @param verbose Whether to print progress messages.
##' @param clusterApply The cluster apply function that is used to distribute the
##'   parallel passes to CPUs in a compute cluster.
##' @param clusterExport A function that is used to export R variables to the nodes of
##'   a CPU cluster, defaults to \code{sfExport}.
##' @return A genetic programming result object that contains a GP population in the
##'   field \code{population}, as well as metadata describing the run parameters.
##'
##' @seealso \code{\link{geneticProgramming}}, \code{\link{summary.geneticProgrammingResult}}, \code{\link{symbolicRegression}}
##' @export
multiNicheGeneticProgramming <- function(fitnessFunction,
                                         stopCondition = makeTimeStopCondition(25),
                                         passStopCondition = makeTimeStopCondition(5),
                                         numberOfNiches = 2,
                                         clusterFunction = groupListConsecutive,
                                         joinFunction = function(niches) Reduce(c, niches),
                                         population = NULL,
                                         populationSize = 100,
                                         eliteSize = ceiling(0.1 * populationSize),
                                         elite = list(),
                                         functionSet = mathFunctionSet,
                                         inputVariables = inputVariableSet("x"),
                                         constantSet = numericConstantSet,
                                         crossoverFunction = crossover,
                                         mutationFunction = NULL,
                                         restartCondition = makeEmptyRestartCondition(),
                                         restartStrategy = makeLocalRestartStrategy(),
                                         searchHeuristic = makeAgeFitnessComplexityParetoGpSearchHeuristic(),
                                         progressMonitor = NULL,
                                         verbose = TRUE,
                                         clusterApply = sfClusterApplyLB,
                                         clusterExport = sfExport) {
  ## Provide default parameters and initialize GP run...
  logmsg <- function(msg, ...) {
    if (verbose)
      message(sprintf(msg, ...))
  }
  quietProgmon <- function(pop, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) NULL
  environment(quietProgmon) <- globalenv() # prevent a possible memory-leak
  progmon <-
    if (verbose) {
      function(pop, objectiveVectors, fitnessFunction, evaluationNumber, stepNumber, bestFitness, timeElapsed, ...) {
        if (!is.null(progressMonitor))
          progressMonitor(pop, objectiveVectors, fitnessFunction, evaluationNumber, stepNumber, bestFitness, timeElapsed, ...)
        if (stepNumber %% 100 == 0)
          logmsg("evolution step %i, fitness evaluations: %i, best fitness: %f, time elapsed: %s",
                 stepNumber, evaluationNumber, bestFitness, formatSeconds(timeElapsed))
      }
    } else if (is.null(progressMonitor)) {
      function(pop, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed, ...) NULL # verbose == FALSE, do not show progress
    } else
      progressMonitor
  mutatefunc <-
    if (is.null(mutationFunction)) {
      function(ind) mutateSubtree(mutateNumericConst(ind),
                                  functionSet, inputVariables, constantSet, mutatesubtreeprob = 0.01)
    } else
      mutationFunction
  pop <-
    if (is.null(population))
      makePopulation(populationSize, functionSet, inputVariables, constantSet)
    else
      population
  popClass <- class(pop)
  variablesToExportToClusterNodes <- c("quietProgmon", "mutatefunc", "restartCondition", "restartStrategy", "crossoverFunction",
                                       "searchHeuristic", "constantSet", "inputVariables", "functionSet", "eliteSize",
                                       "passStopCondition", "fitnessFunction")
  stepNumber <- 1
  evaluationNumber <- 0
  startTime <- proc.time()["elapsed"]
  timeElapsed <- 0
  bestFitness <- Inf

  # Distribute multi-niche GP run to compute cluster...
  clusterExport(list = variablesToExportToClusterNodes)
  passWorker <- function(niche)
    geneticProgramming(fitnessFunction,
                       stopCondition = passStopCondition,
                       population = niche,
                       eliteSize = eliteSize,
                       elite = list(), # each pass starts with an empty elite
                       functionSet = functionSet,
                       inputVariables = inputVariables,
                       constantSet = constantSet,
                       crossoverFunction = crossoverFunction,
                       mutationFunction = mutatefunc,
                       restartCondition = restartCondition,
                       restartStrategy = restartStrategy,
                       searchHeuristic = searchHeuristic,
                       progressMonitor = quietProgmon,
                       verbose = FALSE)
  environment(passWorker) <- globalenv()

  ## Execute multi-niche GP run...
  logmsg("STARTING multi-niche genetic programming evolution run...")
  while (!stopCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                        evaluationNumber = evaluationNumber, timeElapsed = timeElapsed)) {
    logmsg("clustering population into %i niches", numberOfNiches)
    niches <- clusterFunction(pop, numberOfNiches) # cluster population into niches
    for (i in 1:length(niches)) class(niches[[i]]) <- popClass # niches should be of class "gp population"
    logmsg("multi-niche pass with %i niches, evolution steps %i, fitness evaluations: %i, best fitness: %f, time elapsed: %s",
           length(niches), stepNumber, evaluationNumber, bestFitness, formatSeconds(timeElapsed))
    passResults <- clusterApply(niches, passWorker)
    bestFitness <- min(c(as.numeric(Map(function(passResult) passResult$bestFitness, passResults)), bestFitness))
    timeElapsed <- proc.time()["elapsed"] - startTime
    stepNumber <- stepNumber + Reduce(`+`, Map(function(passResult) passResult$stepNumber, passResults))
    evaluationNumber <- evaluationNumber + Reduce(`+`, Map(function(passResult) passResult$evaluationNumber, passResults))
    passResultNiches <- Map(function(passResult) passResult$population, passResults)
    passResultElites <- Map(function(passResult) passResult$elite, passResults)
    allElites <- c(elite, passResultElites)
    pop <- joinFunction(passResultNiches) # join niches again after each pass
    elite <- Reduce(function(e1, e2) joinElites(e1, e2, eliteSize, fitnessFunction), allElites) # join all elites
    class(pop) <- popClass # pop should be of class "gp population"
  }
  logmsg("Multi-niche genetic programming evolution run FINISHED after %i evolution steps, %i fitness evaluations and %s.",
         stepNumber, evaluationNumber, formatSeconds(timeElapsed))

  ## Return GP run result...
  structure(list(fitnessFunction = fitnessFunction,
                 stopCondition = stopCondition,
                 timeElapsed = timeElapsed,
                 stepNumber = stepNumber,
                 evaluationNumber = evaluationNumber,
                 bestFitness = bestFitness,
                 population = pop,
                 elite = elite,
                 functionSet = functionSet,
                 constantSet = constantSet,
                 crossoverFunction = crossoverFunction,
                 mutationFunction = mutatefunc,
                 restartCondition = restartCondition,
                 restartStrategy = restartStrategy), class = "geneticProgrammingResult")
}

##' Symbolic regression via multi-niche standard genetic programming
##'
##' Perform symbolic regression via untyped multi-niche genetic programming.
##' The regression task is specified as a \code{\link{formula}}. Only simple
##' formulas without interactions are supported. The result of the symbolic
##' regression run is a symbolic regression model containing an untyped GP
##' population of model functions.
##'
##' @param formula A \code{\link{formula}} describing the regression task. Only
##'   simple formulas of the form \code{response ~ variable1 + ... + variableN}
##'   are supported at this point in time.
##' @param data A \code{\link{data.frame}} containing training data for the
##'   symbolic regression run. The variables in \code{formula} must match
##'   column names in this data frame.
##' @param stopCondition The stop condition for the evolution main loop. See
##'   \code{makeStepsStopCondition} for details.
##' @param passStopCondition The stop condition for each parallel pass. See
##'   \code{makeStepsStopCondition} for details.
##' @param numberOfNiches The number of niches to cluster the population into.
##' @param clusterFunction The function used to cluster the population into
##'   niches. The first parameter of this function is a GP population, the
##'   second paramater an integer representing the number of niches. Defaults
##'   to \code{\link{groupListConsecutive}}.
##' @param joinFunction The function used to join all niches into a population
##'   again after a round of parallel passes. Defaults to a function that
##'   simply concatenates all niches.
##' @param population The GP population to start the run with. If this parameter
##'   is missing, a new GP population of size \code{populationSize} is created
##'   through random growth.
##' @param populationSize The number of individuals if a population is to be
##'   created.
##' @param eliteSize The number of "elite" individuals to keep. Defaults to
##'  \code{ceiling(0.1 * populationSize)}.
##' @param elite The elite list, must be alist of individuals sorted in ascending
##'   order by their first fitness component.
##' @param individualSizeLimit Individuals with a number of tree nodes that
##'   exceeds this size limit will get a fitness of \code{Inf}.
##' @param penalizeGenotypeConstantIndividuals Individuals that do not contain
##'   any input variables will get a fitness of \code{Inf}.
##' @param functionSet The function set.
##' @param constantSet The set of constant factory functions.
##' @param selectionFunction The selection function to use. Defaults to
##'   tournament selection. See \link{makeTournamentSelection} for details.
##' @param crossoverFunction The crossover function.
##' @param mutationFunction The mutation function.
##' @param restartCondition The restart condition for the evolution main loop. See
##'   \link{makeFitnessStagnationRestartCondition} for details.
##' @param restartStrategy The strategy for doing restarts. See
##'   \link{makeLocalRestartStrategy} for details.
##' @param progressMonitor A function of signature
##'   \code{function(population, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber,
##'   bestFitness, timeElapsed, ...)} to be called with each evolution step. Seach heuristics
##'   may pass additional information via the \code{...} parameter.
##' @param verbose Whether to print progress messages.
##' @param clusterApply The cluster apply function that is used to distribute the
##'   parallel passes to CPUs in a compute cluster.
##' @param clusterExport A function that is used to export R variables to the nodes of
##'   a CPU cluster, defaults to snowfall's \code{sfExport}.
##' @return An symbolic regression model that contains an untyped GP population.
##'
##' @seealso \code{\link{predict.symbolicRegressionModel}}, \code{\link{geneticProgramming}}
##' @export
multiNicheSymbolicRegression <- function(formula, data,
                                         stopCondition = makeTimeStopCondition(25),
                                         passStopCondition = makeTimeStopCondition(5),
                                         numberOfNiches = 2,
                                         clusterFunction = groupListConsecutive,
                                         joinFunction = function(niches) Reduce(c, niches),
                                         population = NULL,
                                         populationSize = 100,
                                         eliteSize = ceiling(0.1 * populationSize),
                                         elite = list(),
                                         individualSizeLimit = 64,
                                         penalizeGenotypeConstantIndividuals = FALSE,
                                         functionSet = mathFunctionSet,
                                         constantSet = numericConstantSet,
                                         selectionFunction = makeTournamentSelection(),
                                         crossoverFunction = crossover,
                                         mutationFunction = NULL,
                                         restartCondition = makeEmptyRestartCondition(),
                                         restartStrategy = makeLocalRestartStrategy(),
                                         progressMonitor = NULL,
                                         verbose = TRUE,
                                         clusterApply = sfClusterApplyLB,
                                         clusterExport = sfExport) {
  ## Match variables in formula to those in data or parent.frame() and
  ## return them in a new data frame. This also expands any '.'
  ## arguments in the formula.  
  mf <- model.frame(formula, data)
  ## Extract list of terms (rhs of ~) in expanded formula
  variableNames <- attr(terms(formula(mf)), "term.labels")
  ## Create inputVariableSet
  inVarSet <- inputVariableSet(list=as.list(variableNames))
  fitFunc <- makeRegressionFitnessFunction(formula(mf), mf, errorMeasure = rmse,
                                           penalizeGenotypeConstantIndividuals = penalizeGenotypeConstantIndividuals,
                                           indsizelimit = individualSizeLimit)
  gpModel <- multiNicheGeneticProgramming(fitFunc, stopCondition, passStopCondition,
                                          numberOfNiches, clusterFunction, joinFunction,
                                          population, populationSize, eliteSize, elite,
                                          functionSet, inVarSet, constantSet, selectionFunction,
                                          crossoverFunction, mutationFunction,
                                          restartCondition, restartStrategy,
                                          progressMonitor, verbose,
                                          clusterApply, clusterExport)
  
  structure(append(gpModel, list(formula = formula(mf))),
                   class = c("symbolicRegressionModel", "geneticProgrammingResult"))
}

##' Clustering Populations for Niching
##'
##' These functions create \code{clusterFunction}s for
##' \code{\link{multiNicheGeneticProgramming}} and \code{\link{multiNicheSymbolicRegression}}.
##' \code{makeHierarchicalClusterFunction} returns a clustering function that uses Ward's
##' agglomerative hierarchical clustering algorithm \code{\link{hclust}}.
##'
##' @param distanceMeasure A distance measure, used for calculating distances between individuals
##'   in a population.
##' @param minNicheSize The minimum number of individuals in each niche.
##' @return A \code{clusterFunction} for clustering populations.
##'
##' @rdname populationClustering
##' @seealso \code{\link{multiNicheGeneticProgramming}}, \code{\link{multiNicheSymbolicRegression}}
##' @export
makeHierarchicalClusterFunction <- function(distanceMeasure = NULL,
                                            minNicheSize = 1) {
  distanceMeasure <- if (is.null(distanceMeasure)) normInducedFunctionDistance(exprVisitationLength) else distanceMeasure
  function(p, numberOfClusters) {
    dm <- customDist(p, distanceMeasure)
    groups <- cutree(hclust(dm, method = "ward"), numberOfClusters)
    Filter(function(niche) length(niche) > minNicheSize, splitList(p, groups))
  }
}

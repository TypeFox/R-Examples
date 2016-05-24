## evolution.R
##   - Functions defining typical evolution main loops,
##     some typical GP function- and constant sets
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Standard typed and untyped genetic programming
##'
##' Perform a standard genetic programming (GP) run. Use \code{geneticProgramming} for
##' untyped genetic programming or \code{typedGeneticProgramming} for typed genetic
##' programming runs. The required argument \code{fitnessFunction} must be supplied with
##' an objective function that assigns a numerical fitness value to an R function. Fitness
##' values are minimized, i.e. smaller values denote higher/better fitness. If a
##' multi-objective \code{selectionFunction} is used, \code{fitnessFunction} return a
##' numerical vector of fitness values. The result of the GP run is a GP result object
##' containing a GP population of R functions. \code{summary.geneticProgrammingResult} can
##' be used to create summary views of a GP result object. During the run, restarts are
##' triggered by the \code{restartCondition}. When a restart is triggered, the restartStrategy
##' is executed, which returns a new population to replace the current one as well as a list of
##' elite individuals. These are added to the runs elite list, where fitter individuals replace
##' individuals with lesser fittness. The runs elite list is always sorted by fitness in
##' ascending order. Only the first component of a multi-criterial fitness counts in this
##' sorting. After a GP run, the population is inserted into the elite list. The elite list
##' is returned as part of the GP result object.
##'
##' @param fitnessFunction In case of a single-objective selection function,
##'   \code{fitnessFunction} must be a single function that assigns a
##'   numerical fitness value to a GP individual represented as a R function.
##'   Smaller fitness values mean higher/better fitness. If a multi-objective
##'   selection function is used, \code{fitnessFunction} must return a numerical
##'   vector of fitness values.
##' @param type The range type of the individual functions. This parameter
##'   only applies to \code{typedGeneticProgramming}.
##' @param stopCondition The stop condition for the evolution main loop. See
##'   code \code{makeStepsStopCondition} for details.
##' @param population The GP population to start the run with. If this parameter
##'   is missing, a new GP population of size \code{populationSize} is created
##'   through random growth.
##' @param populationSize The number of individuals if a population is to be
##'   created.
##' @param eliteSize The number of elite individuals to keep. Defaults to
##'  \code{ceiling(0.1 * populationSize)}.
##' @param elite The elite list, must be alist of individuals sorted in ascending
##'   order by their first fitness component.
##' @param functionSet The function set.
##' @param inputVariables The input variable set.
##' @param constantSet The set of constant factory functions.
##' @param crossoverFunction The crossover function.
##' @param mutationFunction The mutation function.
##' @param restartCondition The restart condition for the evolution main loop. See
##'   \link{makeEmptyRestartCondition} for details.
##' @param restartStrategy The strategy for doing restarts. See
##'   \link{makeLocalRestartStrategy} for details.
##' @param searchHeuristic The search-heuristic (i.e. optimization algorithm) to use
##'   in the search of solutions. See the documentation for \code{searchHeuristics} for
##'   available algorithms.
##' @param breedingFitness A "breeding" function. This function is applied after
##'   every stochastic operation \emph{Op} that creates or modifies an individal
##'   (typically, \emph{Op} is a initialization, mutation, or crossover operation). If
##'   the breeding function returns \code{TRUE} on the given individual, \emph{Op} is
##'   considered a success. If the breeding function returns \code{FALSE}, \emph{Op}
##'   is retried a maximum of \code{breedingTries} times. If this maximum number of
##'   retries is exceeded, the result of the last try is considered as the result of
##'   \emph{Op}. In the case the breeding function returns a numeric value, the breeding
##'   is repeated \code{breedingTries} times and the individual with the lowest breeding
##'   fitness is considered the result of \emph{Op}.
##' @param breedingTries In case of a boolean \code{breedingFitness} function, the
##'   maximum number of retries. In case of a numerical \code{breedingFitness} function,
##'   the number of breeding steps. Also see the documentation for the \code{breedingFitness}
##'   parameter. Defaults to \code{50}.
##' @param extinctionPrevention When set to \code{TRUE}, the initialization and
##'   selection steps will try to prevent duplicate individuals
##'   from occurring in the population. Defaults to \code{FALSE}, as this
##'   operation might be expensive with larger population sizes.
##' @param archive If set to \code{TRUE}, all GP individuals evaluated are stored in an
##'   archive list \code{archiveList} that is returned as part of the result of this function. 
##' @param progressMonitor A function of signature
##'   \code{function(population, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber,
##'   bestFitness, timeElapsed, ...)} to be called with each evolution step. Seach heuristics
##'   may pass additional information via the \code{...} parameter.
##' @param verbose Whether to print progress messages.
##' @return A genetic programming result object that contains a GP population in the
##'   field \code{population}, as well as metadata describing the run parameters.
##'
##' @seealso \code{\link{summary.geneticProgrammingResult}}, \code{\link{symbolicRegression}}
##' @rdname geneticProgramming
##' @export
geneticProgramming <- function(fitnessFunction,
                               stopCondition = makeTimeStopCondition(5),
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
                               searchHeuristic = makeAgeFitnessComplexityParetoGpSearchHeuristic(lambda = ceiling(0.5 * populationSize)),
                               breedingFitness = function(individual) TRUE,
                               breedingTries = 50,
                               extinctionPrevention = FALSE,
                               archive = FALSE,
                               progressMonitor = NULL,
                               verbose = TRUE) {
  ## Provide default parameters and initialize GP run...
  logmsg <- function(msg, ...) {
    if (verbose)
      message(sprintf(msg, ...))
  }
  nonVerboseProgmon <- if (is.null(progressMonitor)) {
    function(pop, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed, ...) NULL
  } else {
    progressMonitor
  }
  progmon <- if (verbose) {
    function(pop, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed, ...) {
      nonVerboseProgmon(pop, objectiveVectors, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed, ...)
      if (stepNumber %% 100 == 0) {
        logmsg("evolution step %i, fitness evaluations: %i, best fitness: %f, time elapsed: %s",
               stepNumber, evaluationNumber, bestFitness, formatSeconds(timeElapsed))
      }
    }
  } else {
    nonVerboseProgmon
  }
  mutatefunc <-
    if (is.null(mutationFunction)) {
      function(ind) mutateSubtree(mutateNumericConst(ind),
                                  functionSet, inputVariables, constantSet, mutatesubtreeprob = 0.1,
                                  breedingFitness = breedingFitness, breedingTries = breedingTries)
    } else
      mutationFunction
  pop <-
    if (is.null(population))
      makePopulation(populationSize, functionSet, inputVariables, constantSet,
                     extinctionPrevention = extinctionPrevention,
                     breedingFitness = breedingFitness, breedingTries = breedingTries)
    else
      population

  ## Execute search-heuristic...
  result <- searchHeuristic(logFunction = logmsg, stopCondition = stopCondition,
                            pop = pop, fitnessFunction = fitnessFunction,
                            mutationFunction = mutatefunc, crossoverFunction = crossoverFunction,
                            functionSet = functionSet, inputVariables = inputVariables, constantSet = constantSet,
                            archive = archive, extinctionPrevention = extinctionPrevention,
                            elite = elite, eliteSize = eliteSize,
                            restartCondition = restartCondition, restartStrategy = restartStrategy,
                            breedingFitness = breedingFitness, breedingTries = breedingTries,
                            progressMonitor = progmon)

  ## Return GP run result...
  structure(list(fitnessFunction = fitnessFunction,
                 stopCondition = stopCondition,
                 timeElapsed = result$timeElapsed,
                 stepNumber = result$stepNumber,
                 evaluationNumber = result$evaluationNumber,
                 bestFitness = result$bestFitness,
                 population = result$population,
                 fitnessValues = result$fitnessValues,
                 elite = result$elite,
                 functionSet = functionSet,
                 constantSet = constantSet,
                 crossoverFunction = crossoverFunction,
                 mutationFunction = mutatefunc,
                 restartCondition = restartCondition,
                 breedingFitness = breedingFitness,
                 breedingTries = breedingTries,
                 extinctionPrevention = extinctionPrevention,
                 archive = archive,
                 archiveList = result$archiveList,
                 restartStrategy = restartStrategy,
                 searchHeuristicResults = result$searchHeuristicResults),
            class = "geneticProgrammingResult")
}

##' @rdname geneticProgramming
##' @export
typedGeneticProgramming <- function(fitnessFunction,
                                    type,
                                    stopCondition = makeTimeStopCondition(5),
                                    population = NULL,
                                    populationSize = 100,
                                    eliteSize = ceiling(0.1 * populationSize),
                                    elite = list(),
                                    functionSet,
                                    inputVariables,
                                    constantSet,
                                    crossoverFunction = crossoverTyped,
                                    mutationFunction = NULL,
                                    restartCondition = makeEmptyRestartCondition(),
                                    restartStrategy = makeLocalRestartStrategy(populationType = type),
                                    searchHeuristic = makeAgeFitnessComplexityParetoGpSearchHeuristic(),
                                    breedingFitness = function(individual) TRUE,
                                    breedingTries = 50,
                                    extinctionPrevention = FALSE,
                                    archive = FALSE,
                                    progressMonitor = NULL,
                                    verbose = TRUE) {
  if (is.null(type)) stop("typedGeneticProgramming: Type must not be NULL.")
  pop <-
    if (is.null(population))
      makeTypedPopulation(populationSize, type, functionSet, inputVariables, constantSet,
                          extinctionPrevention = extinctionPrevention,
                          breedingFitness = breedingFitness, breedingTries = breedingTries)
    else
      population
  mutatefunc <-
    if (is.null(mutationFunction)) {
      function(ind) mutateSubtreeTyped(mutateFuncTyped(mutateNumericConstTyped(ind),
                                                       functionSet, mutatefuncprob = 0.01),
                                       functionSet, inputVariables, constantSet,
                                       mutatesubtreeprob = 0.01,
                                       breedingFitness = breedingFitness, breedingTries = breedingTries)
    } else mutationFunction
  geneticProgramming(fitnessFunction, stopCondition = stopCondition, population = pop,
                     populationSize = populationSize, eliteSize = eliteSize, elite = elite,
                     functionSet = functionSet,
                     inputVariables = inputVariables,
                     constantSet = constantSet,
                     crossoverFunction = crossoverFunction, mutationFunction = mutatefunc,
                     restartCondition = restartCondition, restartStrategy = restartStrategy,
                     searchHeuristic = searchHeuristic,
                     breedingFitness = breedingFitness, breedingTries = breedingTries,
                     extinctionPrevention = extinctionPrevention,
                     archive = archive,
                     progressMonitor = progressMonitor, verbose = verbose)
}

##' Join elite lists
##'
##' Inserts a list of new individuals into an elite list, replacing the worst individuals
##' in this list to make place, if needed.
##'
##' @param individuals The list of individuals to insert.
##' @param elite The list of elite individuals to insert \code{individuals} into. This
##'   list must be sorted by fitness in ascending order, i.e. lower fitnesses first.
##' @param eliteSize The maximum size of the \code{elite}.
##' @param fitnessFunction The fitness function.
##' @return The \code{elite} with \code{individuals} inserted, sorted by fitness in
##'   ascending order, i.e. lower fitnesses first.
joinElites <- function(individuals, elite, eliteSize, fitnessFunction) {
  allIndividuals <- c(individuals, elite)
  allFitnesses <- as.numeric(Map(function(ind) fitnessFunction(ind)[1], allIndividuals))
  allIndividualsSorted <- allIndividuals[order(allFitnesses, decreasing = FALSE)]
  newElite <- if (length(allIndividualsSorted) > eliteSize)
    allIndividualsSorted[1:eliteSize]
  else
    allIndividualsSorted
  newElite
}

##' Summary reports of genetic programming run result objects
##'
##' Create a summary report of a genetic programming result object as returned by
##' \code{\link{geneticProgramming}} or \code{\link{symbolicRegression}}, for
##' example.
##'
##' @param object The genetic programming run result object to report on.
##' @param reportFitness Whether to report detailed fitness values of each individual
##'   in the result population. Note that calculating fitness values may take
##'   a long time. Defaults to \code{FALSE}. Either way, basic fitness
##'   values for each individual is reported.
##' @param orderByFitness Whether the report of the result population should be
##'   ordered by fitness. This does not have an effect if \code{reportFitness}
##'   is set to \code{FALSE}. Defaults to \code{TRUE}.
##' @param ... Ignored in this summary function.
##'
##' @seealso \code{\link{geneticProgramming}}, \code{\link{symbolicRegression}}
##' @method summary geneticProgrammingResult
##' @S3method summary geneticProgrammingResult
##' @export
summary.geneticProgrammingResult <- function(object, reportFitness = FALSE, orderByFitness = TRUE, ...) {
  reportPopulation <- function(individualFunctions, individualFitnessValues) {
    individualFunctionsAsStrings <- Map(function(f) Reduce(function(a, b) paste(a, b, sep=""),
                                                           deparse(f)),
                                        individualFunctions)
    report <- cbind(1:length(individualFunctions), individualFunctions, individualFunctionsAsStrings, individualFitnessValues)
    colnames(report) <- c("Individual Index", "Individual Function", "(as String)", "Fitness Value")
    if (reportFitness) {
      fitnessList <- lapply(individualFunctions, object$fitnessFunction)
      fitnessesLength <- length(fitnessList)
      fitnessDimemsion <- length(fitnessList[[1]])
      flatFitnesses <- Reduce(c, fitnessList)
      fitnessMatrix <- matrix(as.list(flatFitnesses), ncol = fitnessDimemsion)
      if (is.null(colnames(fitnessMatrix))) colnames(fitnessMatrix) <- paste("Fitness", 1:fitnessDimemsion)
      report <- cbind(report, fitnessMatrix)
      if (orderByFitness) {
        rawFitnessMatrix <- matrix(flatFitnesses, ncol = fitnessDimemsion)
        report <- report[do.call(order, split(rawFitnessMatrix, col(rawFitnessMatrix))),]
      }
    }
    report
  }
  list(population = reportPopulation(object$population, object$fitnessValues), elite = reportPopulation(object$elite))
}

##' Evolution restart conditions
##'
##' Evolution restart conditions are predicates (functions that return a single logical value)
##' of the signature \code{function(population, fitnessFunction, stepNumber, evaluationNumber,
##' bestFitness, timeElapsed)}. They are used to decide when to restart a GP evolution run that
##' might be stuck in a local optimum. Evolution restart conditions are objects of the same
##' type and class as evolution stop conditions. They may be freely substituted for each other.
##'
##' \code{makeEmptyRestartCondition} creates a restart condition that is never fulfilled, i.e.
##' restarts will never occur.
##' \code{makeStepLimitRestartCondition} creates a restart condition that holds if the
##' number if evolution steps is an integer multiple of a given step limit.
##' restarts will never occur.
##' \code{makeFitnessStagnationRestartCondition} creates a restart strategy that holds if the
##' standard deviation of a last \code{fitnessHistorySize} best fitness values falls below
##' a given \code{fitnessStandardDeviationLimit}.
##' \code{makeFitnessDistributionRestartCondition} creates a restart strategy that holds
##' if the standard deviation of the fitness values of the individuals in the current
##' population falls below a given \code{fitnessStandardDeviationLimit}.
##'
##' @param stepLimit The step limit for \code{makeStepLimitRestartCondition}.
##' @param fitnessHistorySize The number of past best fitness values to look at when calculating
##'   the best fitness standard deviation for \code{makeFitnessStagnationRestartCondition}.
##' @param testFrequency The frequency to test for the restart condition, in evolution
##'   steps. This parameter is mainly used with restart condititions that are expensive to
##'   calculate.
##' @param fitnessStandardDeviationLimit The best fitness standard deviation limit for
##'   \code{makeFitnessStagnationRestartCondition}.
##'
##' @rdname evolutionRestartConditions
##' @export
makeEmptyRestartCondition <- function() {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) FALSE
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionRestartConditions
##' @export
makeStepLimitRestartCondition <- function(stepLimit = 10) {
  function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) {
    stepNumber %% stepLimit == 0
  }
}

##' @rdname evolutionRestartConditions
##' @export
makeFitnessStagnationRestartCondition <- function(fitnessHistorySize = 100,
                                                  testFrequency = 10,
                                                  fitnessStandardDeviationLimit = 0.000001) {
  history <- rep(Inf, fitnessHistorySize)
  # TODO this restart condition is broken by design, because best fitness is not affected by restarts!
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) {
    if (stepNumber %% testFrequency == 0) {
      history <<- c(history[2:fitnessHistorySize], bestFitness)
      historySd <- sd(history)
      !is.nan(historySd) && historySd <= fitnessStandardDeviationLimit # restart if sd drop below limit
    } else FALSE
  }
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionRestartConditions
##' @export
makeFitnessDistributionRestartCondition <- function(testFrequency = 100,
                                                    fitnessStandardDeviationLimit = 0.000001) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) {
    if (stepNumber %% testFrequency == 0) {
      fitnessValues <- as.numeric(Map(function(ind) fitnessFunction(ind)[1], pop))
      fitnessSd <- sd(fitnessValues)
      !is.nan(fitnessSd) && fitnessSd <= fitnessStandardDeviationLimit # restart if sd drop below limit
    } else FALSE
  }
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' Evolution restart strategies
##'
##' Evolution restart strategies are functions of the signature \code{function(fitnessFunction,
##' population, populationSize, functionSet, inputVariables, constantSet)} that return a list of
##' two obtjects: First, a population that replace the run's current population. Second, a list
##' of elite individuals to keep.
##'
##' \code{makeLocalRestartStrategy} creates a restart strategy that replaces all individuals with
##' new individuals. The single best individual is returned as the elite. When using a
##' multi-criterial fitness function, only the first component counts in the fitness sorting.
##'
##' @param populationType The sType of the replacement individuals, defaults to \code{NULL} for
##'   creating untyped populations.
##' @param extinctionPrevention Whether to surpress duplicate individuals in newly initialized
##'   populations. See \code{\link{geneticProgramming}} for details.
##' @param breedingFitness A breeding function. See the documentation for
##'   \code{\link{geneticProgramming}} for details.
##' @param breedingTries The number of breeding steps.
##'
##' @rdname evolutionRestartStrategies
##' @export
makeLocalRestartStrategy <- function(populationType = NULL,
                                     extinctionPrevention = FALSE,
                                     breedingFitness = function(individual) TRUE,
                                     breedingTries = 50) {
  restartStrategy <- function(fitnessFunction, population, populationSize, functionSet, inputVariables, constantSet) {
    fitnessValues <- as.numeric(Map(function(ind) fitnessFunction(ind)[1], population))
    bestInd <- population[[which.min(fitnessValues)]]
    restartedPopulation <-
      if (is.null(populationType))
        makePopulation(populationSize, functionSet, inputVariables, constantSet,
                       extinctionPrevention = extinctionPrevention,
                       breedingFitness = breedingFitness, breedingTries = breedingTries)
      else
        makeTypedPopulation(populationSize, populationType, functionSet, inputVariables, constantSet,
                            extinctionPrevention = extinctionPrevention,
                            breedingFitness = breedingFitness, breedingTries = breedingTries)
    list(population = restartedPopulation, elite = list(bestInd))
  }
  class(restartStrategy) <- c("restartStrategy", "function")
  restartStrategy
}

##' Evolution stop conditions
##'
##' Evolution stop conditions are predicates (functions that return a single logical value)
##' of the signature \code{function(population, stepNumber, evaluationNumber, bestFitness,
##' timeElapsed)}.
##' They are used to decide when to finish a GP evolution run. Stop conditions must be members
##' of the S3 class \code{c("stopCondition", "function")}. They can be combined using the
##' functions \code{andStopCondition}, \code{orStopCondition} and \code{notStopCondition}.
##'
##' \code{makeStepsStopCondition} creates a stop condition that is fulfilled if the number
##' of evolution steps exceeds a given limit.
##' \code{makeEvaluationsStopCondition} creates a stop condition that is fulfilled if the
##' number of fitness function evaluations exceeds a given limit.
##' \code{makeFitnessStopCondition} creates a stop condition that is fulfilled if the
##' number best fitness seen in an evaluation run undercuts a certain limit.
##' \code{makeTimeStopCondition} creates a stop condition that is fulfilled if the run time
##' (in seconds) of an evolution run exceeds a given limit.
##'
##' @param stepLimit The maximum number of evolution steps for \code{makeStepsStopCondition}.
##' @param evaluationLimit The maximum number of fitness function evaluations for
##'   \code{makeEvaluationsStopCondition}.
##' @param fitnessLimit The minimum fitness for \code{makeFitnessStopCondition}.
##' @param timeLimit The maximum runtime in seconds for \code{makeTimeStopCondition}.
##' @param e1 A stop condition.
##' @param e2 A stop condition.
##'
##' @rdname evolutionStopConditions
##' @export
makeStepsStopCondition <- function(stepLimit) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) stepNumber >= stepLimit
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionStopConditions
##' @export
makeEvaluationsStopCondition <- function(evaluationLimit) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) evaluationNumber >= evaluationLimit
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionStopConditions
##' @export
makeFitnessStopCondition <- function(fitnessLimit) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) bestFitness <= fitnessLimit
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionStopConditions
##' @export
makeTimeStopCondition <- function(timeLimit) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) timeElapsed >= timeLimit
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionStopConditions
##' @export
andStopCondition <- function(e1, e2) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
    e1(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) && e2(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionStopConditions
##' @export
orStopCondition <- function(e1, e2) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
    e1(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed) || e2(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' @rdname evolutionStopConditions
##' @export
notStopCondition <- function(e1) {
  stopCondition <- function(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
    !e1(pop, fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
  class(stopCondition) <- c("stopCondition", "function")
  stopCondition
}

##' Some simple arithmetic and logic functions for use in GP expressions
##'
##' \code{safeDivide} a division operator that returns 0 if the divisor is 0.
##' \code{safeLn} a natural logarithm operator that return 0 if its argument is less
##'   then 0.
##' \code{ln} is the natural logarithm.
##' \code{positive} returns true if its argument is greater then 0.
##' \code{ifPositive} returns its second argument if its first argument is positive,
##'   otherwise its third argument.
##' \code{ifThenElse} returns its second argument if its first argument is \code{TRUE},
##'   otherwise its third argument.
##'
##' @param a A numeric value.
##' @param b A numeric value.
##' @param x A numeric value.
##' @param thenbranch The element to return when \code{x} is \code{TRUE}.
##' @param elsebranch The element to return when \code{x} is \code{FALSE}.
##'
##' @rdname safeGPfunctions
##' @export
safeDivide <- function(a, b) ifelse(b == 0, b, a / b)

##' @rdname safeGPfunctions
##' @export
safeSqroot <- function(a) sqrt(ifelse(a < 0, 0, a))

##' @rdname safeGPfunctions
##' @export
safeLn <- function(a) ifelse(a < 0, 0, log(a))

##' @rdname safeGPfunctions
##' @export
ln <- function(a) log(a)

##' @rdname safeGPfunctions
##' @export
positive <- function(x) x > 0

##' @rdname safeGPfunctions
##' @export
ifPositive <- function(x, thenbranch, elsebranch) ifelse(x > 0, thenbranch, elsebranch)

##' @rdname safeGPfunctions
##' @export
ifThenElse <- function(x, thenbranch, elsebranch) ifelse(x, thenbranch, elsebranch)

# TODO hack to make roxygen2 work with external function references...
functionSet <- function(...) NULL
constantFactorySet <- function(...) NULL
# ... .

##' Default function- and constant factory sets for Genetic Programming
##'
##' \code{arithmeticFunctionSet} is an untyped function set containing the functions
##' "+", "-", "*", and "/".
##' \code{expLogFunctionSet} is an untyped function set containing the functions
##' "sqrt", "exp", and "ln".
##' \code{trigonometricFunctionSet} is an untyped function set containing the functions
##' "sin", "cos", and "tan".
##' \code{mathFunctionSet} is an untyped function set containing all of the above functions.
##'
##' \code{numericConstantSet} is an untyped constant factory set containing a single
##' constant factory that creates numeric constants via calls to \code{runif(1, -1, 1)}.
##'
##' Note that these objects are initialized in the RGP package's \code{.onAttach} function.
##'
##' @rdname defaultGPFunctionAndConstantSets
##' @export
arithmeticFunctionSet <- NULL

##' @rdname defaultGPFunctionAndConstantSets
##' @export
expLogFunctionSet <- NULL 

##' @rdname defaultGPFunctionAndConstantSets
##' @export
#trigonometricFunctionSet <- functionSet("sin", "cos", "tan")

##' @rdname defaultGPFunctionAndConstantSets
##' @export
#mathFunctionSet <- c(arithmeticFunctionSet, expLogFunctionSet, trigonometricFunctionSet)

##' @rdname defaultGPFunctionAndConstantSets
##' @export
#numericConstantSet <- constantFactorySet(function() runif(1, -1, 1))

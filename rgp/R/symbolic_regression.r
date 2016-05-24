## symbolic_regression.R
##   - Tools for symbolic regression using the R formula interface
##
## RGP - a GP system for R
## 2010-2011 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Symbolic regression via untyped standard genetic programming
##'
##' Perform symbolic regression via untyped genetic programming. The regression
##' task is specified as a \code{\link{formula}}. Only simple formulas without
##' interactions are supported. The result of the symbolic regression run is a
##' symbolic regression model containing an untyped GP population of model
##' functions.
##'
##' @param formula A \code{\link{formula}} describing the regression task. Only
##'   simple formulas of the form \code{response ~ variable1 + ... + variableN}
##'   are supported at this point in time.
##' @param data A \code{\link{data.frame}} containing training data for the
##'   symbolic regression run. The variables in \code{formula} must match
##'   column names in this data frame.
##' @param stopCondition The stop condition for the evolution main loop. See
##'   \link{makeStepsStopCondition} for details.
##' @param population The GP population to start the run with. If this parameter
##'   is missing, a new GP population of size \code{populationSize} is created
##'   through random growth.
##' @param populationSize The number of individuals if a population is to be
##'   created.
##' @param eliteSize The number of elite individuals to keep. Defaults to
##'  \code{ceiling(0.1 * populationSize)}.
##' @param elite The elite list, must be alist of individuals sorted in ascending
##'   order by their first fitness component.
##' @param extinctionPrevention When set to \code{TRUE}, the initialization and
##'   selection steps will try to prevent duplicate individuals
##'   from occurring in the population. Defaults to \code{FALSE}, as this
##'   operation might be expensive with larger population sizes.
##' @param archive If set to \code{TRUE}, all GP individuals evaluated are stored in an
##'   archive list \code{archiveList} that is returned as part of the result of this function. 
##' @param individualSizeLimit Individuals with a number of tree nodes that
##'   exceeds this size limit will get a fitness of \code{Inf}.
##' @param penalizeGenotypeConstantIndividuals Individuals that do not contain
##'   any input variables will get a fitness of \code{Inf}.
##' @param subSamplingShare The share of fitness cases \deqn{s} sampled for
##'   evaluation with each function evaluation. \deqn{0 < s \leq 1} must
##'   hold, defaults to \code{1.0}.
##' @param functionSet The function set.
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
##' @param errorMeasure A function to use as an error measure, defaults to RMSE.
##' @param progressMonitor A function of signature
##'   \code{function(population, fitnessValues, fitnessFunction, stepNumber, evaluationNumber,
##'   bestFitness, timeElapsed)} to be called with each evolution step.
##' @param envir The R environment to evaluate individuals in, defaults to
##'   \code{parent.frame()}.
##' @param verbose Whether to print progress messages.
##' @return An symbolic regression model that contains an untyped GP population.
##'
##' @seealso \code{\link{predict.symbolicRegressionModel}}, \code{\link{geneticProgramming}}
##' @export
symbolicRegression <- function(formula, data,
                               stopCondition = makeTimeStopCondition(5),
                               population = NULL,
                               populationSize = 100,
                               eliteSize = ceiling(0.1 * populationSize),
                               elite = list(),
                               extinctionPrevention = FALSE,
                               archive = FALSE,
                               individualSizeLimit = 64,
                               penalizeGenotypeConstantIndividuals = FALSE,
                               subSamplingShare = 1.0,
                               functionSet = mathFunctionSet,
                               constantSet = numericConstantSet,
                               crossoverFunction = NULL,
                               mutationFunction = NULL,
                               restartCondition = makeEmptyRestartCondition(),
                               restartStrategy = makeLocalRestartStrategy(),
                               searchHeuristic = makeAgeFitnessComplexityParetoGpSearchHeuristic(),
                               breedingFitness = function(individual) TRUE,
                               breedingTries = 50,
                               errorMeasure = rmse,
                               progressMonitor = NULL,
                               envir = parent.frame(),
                               verbose = TRUE) {
  symbolicRegressionModel <- dataDrivenGeneticProgramming(formula, data, makeRegressionFitnessFunction,
                                                          list(envir = envir,
                                                               errorMeasure = errorMeasure,
                                                               penalizeGenotypeConstantIndividuals = penalizeGenotypeConstantIndividuals,
                                                               subSamplingShare = subSamplingShare,
                                                               indsizelimit = individualSizeLimit),
                                                          stopCondition = stopCondition,
                                                          population = population,
                                                          populationSize = populationSize,
                                                          eliteSize = eliteSize,
                                                          elite = elite,
                                                          extinctionPrevention = extinctionPrevention,
                                                          archive = archive,
                                                          functionSet = functionSet,
                                                          constantSet = constantSet,
                                                          crossoverFunction = crossoverFunction,
                                                          mutationFunction = mutationFunction,
                                                          restartCondition = restartCondition,
                                                          restartStrategy = restartStrategy,
                                                          searchHeuristic = searchHeuristic,
                                                          breedingFitness = breedingFitness,
                                                          breedingTries = breedingTries,
                                                          progressMonitor = progressMonitor,
                                                          verbose = verbose)

  class(symbolicRegressionModel) <- c("symbolicRegressionModel", class(symbolicRegressionModel))
  symbolicRegressionModel  
}

##' Predict method for symbolic regression models
##'
##' Predict values via a model function from a population of model functions
##' generated by symbolic regression.
##'
##' @param object A model created by \code{\link{symbolicRegression}}.
##' @param newdata A \code{\link{data.frame}} containing input data for the
##'   symbolic regression model. The variables in \code{object$formula} must match
##'   column names in this data frame.
##' @param model The numeric index of the model function in \code{object$population}
##'   to use for prediction or \code{"BEST"} to use the model function with the best
##'   training fitness.
##' @param detailed Whether to add metadata to the prediction object returned.
##' @param ... Ignored in this \code{predict} method.
##' @return A vector of predicted values or, if \code{detailed} is \code{TRUE}, a
##'   list of the following elements:
##'   \code{model} the model used in this prediction
##'   \code{response} a matrix of predicted versus respone values
##'   \code{RMSE} the RMSE between the real and predicted response
##'
##' @method predict symbolicRegressionModel
##' @S3method predict symbolicRegressionModel
##' @export
predict.symbolicRegressionModel <- function(object, newdata, model = "BEST", detailed = FALSE, ...) {
  ind <- if (model == "BEST") {
    trainingFitnessSortedPopulation <- sortBy(object$population, object$fitnessFunction)
    trainingFitnessSortedPopulation[[1]]
  } else object$population[[model]]
  data <- if (any(is.na(newdata))) {
    dataWithoutNAs <- na.omit(newdata)
    warning(sprintf("removed %i data rows containing NA values", length(attr(dataWithoutNAs, "na.action"))))
    dataWithoutNAs
  } else newdata
  
  formulaVars <- as.list(attr(terms(object$formula), "variables")[-1])
  responseVariable <- formulaVars[[1]]
  explanatoryVariables <- formulaVars[-1]
  trueResponse <- with(data, eval(responseVariable))
  explanatories <- with(data, lapply(explanatoryVariables, eval))
  ysind <- do.call(ind, explanatories) # vectorized evaluation
  errorind <- rmse(trueResponse, ysind)
  
  if (detailed) {
    predictedVersusReal <- cbind(ysind, trueResponse)
    colnames(predictedVersusReal) <- c("predicted", "real")
    list(model = ind, response = predictedVersusReal, RMSE = errorind)
  } else ysind
}

##' Create a fitness function for symbolic regression
##'
##' Creates a fitness function that calculates an error measure with
##' respect to a given set of data variables. A simplified version of
##' the formula syntax is used to describe the regression task. When
##' an \code{indsizelimit} is given, individuals exceeding this limit
##' will receive a fitness of \code{Inf}.
##'
##' @param formula A formula object describing the regression task.
##' @param data An optional data frame containing the variables in the
##'   model.
##' @param envir The R environment to evaluate individuals in.
##' @param errorMeasure A function to use as an error measure, defaults to RMSE.
##' @param indsizelimit Individuals exceeding this size limit will get
##'   a fitness of \code{Inf}.
##' @param penalizeGenotypeConstantIndividuals Individuals that do not
##'   contain any input variables will get a fitness of \code{Inf}.
##' @param subSamplingShare The share of fitness cases \deqn{s} sampled for
##'   evaluation with each function evaluation. \deqn{0 < s \leq 1} must
##'   hold, defaults to \code{1.0}.
##' @return A fitness function to be used in symbolic regression.
##' @export
makeRegressionFitnessFunction <- function(formula, data, envir,
                                          errorMeasure = rmse,
                                          indsizelimit = NA,
                                          penalizeGenotypeConstantIndividuals = FALSE,
                                          subSamplingShare = 1.0) {
  data <- if (any(is.na(data))) {
    dataWithoutNAs <- na.omit(data)
    warning(sprintf("removed %i data rows containing NA values",
                    length(attr(dataWithoutNAs, "na.action"))))
    dataWithoutNAs
  } else data
  formulaVars <- as.list(attr(terms(formula), "variables")[-1])
  responseVariable <- formulaVars[[1]]
  explanatoryVariables <- formulaVars[-1]

  allTrueResponse <- eval(responseVariable, envir = data)
  allExplanatories <- lapply(explanatoryVariables, eval, envir = data)

  numberOfFitnessCases <- length(allTrueResponse)
  sampleIndices <- if (subSamplingShare < 1.0)
    sample(1:numberOfFitnessCases, subSamplingShare * numberOfFitnessCases, replace = FALSE)
  else
    1:numberOfFitnessCases

  explanatories <- Map(function(explanatory) explanatory[sampleIndices], allExplanatories)
  trueResponse <- allTrueResponse[sampleIndices]

  function(ind) {
    ysind <- do.call(ind, explanatories, envir = envir) # vectorized fitness-case evaluation
    errorind <- errorMeasure(trueResponse, ysind)
    if (!is.na(indsizelimit) && funcSize(ind) > indsizelimit)
      Inf # individual size limit exceeded
    else if (is.na(errorind) || is.nan(errorind))
      Inf # error value is NA or NaN
    else if (penalizeGenotypeConstantIndividuals
             && is.empty(inputVariablesOfIndividual(ind, explanatoryVariables)))
      Inf # individual does not contain any input variables
    else errorind
  }
}

##' Variable Presence Maps
##'
##' Counts the number of input variables (formal arguments) present in the
##' body of a individual function. Applied to a population of individuals,
##' this information is useful to identify driving variables in a modelling
##' task.
##' \code{functionVariablePresenceMap} returns a (one row) variable
##' presence map for a function, \code{populationVariablePresenceMap}
##' returns a variable presence map for a population of RGP individuals
##' (a list of R functions).
##'
##' @param f A R function to return a variable presence map for.
##' @param pop A RGP population to return a variable presence map for.
##' @return A data frame with variables (formal parameters) in the columns,
##'   individuals (function) in the rows and variable counts in the cells.
##'
##' @rdname variablePresenceMaps 
##' @export
functionVariablePresenceMap <- function(f) {
  variablePresence <- as.list(formals(f))
  variablePresence <- Map(function(x) 0, variablePresence) # set each element to 0

  MapExpressionLeafs(function(leaf) if (is.name(leaf)) {
    variablePresence[[as.character(leaf)]] <<- variablePresence[[as.character(leaf)]] + 1
  }, body(f))

  as.data.frame(variablePresence)
}

##' @rdname variablePresenceMaps 
##' @export
populationVariablePresenceMap <- function(pop)
  Reduce(rbind, Map(functionVariablePresenceMap, pop))


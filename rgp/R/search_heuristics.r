## search_heuristics.R
##   - Functions defining search heuristics (i.e. algorithmic frameworks) 
##     for Genetic Programming 
##
## RGP - a GP system for R
## 2011 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Tiny GP Search Heuristic for RGP
##'
##' The search-heuristic, i.e. the concrete GP search algorithm, is a modular component of RGP.
##' \code{makeTinyGpSearchHeuristic} creates an RGP search-heuristic that mimics the search heuristic
##' implemented in Riccardo Poli's TinyGP system.
##'
##' @param crossoverProbability The crossover probability for search-heuristics that support
##'   this setting (i.e. TinyGP). Defaults to \code{0.9}.
##' @param tournamentSize The size of TinyGP's selection tournaments.
##' @return An RGP search heuristic.
##'
##' @export
makeTinyGpSearchHeuristic <- function(crossoverProbability = 0.9, tournamentSize = 2)
function(logFunction, stopCondition, pop, fitnessFunction,
         mutationFunction, crossoverFunction,
         functionSet, inputVariables, constantSet,
         archive, extinctionPrevention,
         elite, eliteSize,
         restartCondition, restartStrategy,
         breedingFitness, breedingTries,
         progressMonitor) {
  logFunction("STARTING genetic programming evolution run (TinyGP search-heuristic) ...")
  
  ## Global variables...
  popSize <- length(pop)
  fitnessValues <- sapply(pop, fitnessFunction)

  ## Initialize statistic counters...
  stepNumber <- 1
  evaluationNumber <- 0
  timeElapsed <- 0
  archiveList <- list() # the archive of all individuals selected in this run, only used if archive == TRUE
  archiveIndexOf <- function(archive, individual)
    Position(function(a) identical(body(a$individual), body(individual)), archive)
  bestFitness <- min(fitnessValues) # best fitness value seen in this run, if multi-criterial, only the first component counts
  startTime <- proc.time()["elapsed"]

  ## Execute GP run...
  while (!stopCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                        evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
    child <- NULL
    if (runif(1) < crossoverProbability) {
      # create child via crossover
      motherIndex <- tournament(fitnessValues, popSize, tournamentSize)
      fatherIndex <- tournament(fitnessValues, popSize, tournamentSize)
      child <- crossoverFunction(pop[[motherIndex]], pop[[fatherIndex]],
                                 breedingFitness = breedingFitness,
                                 breedingTries = breedingTries)
    } else {
      # create child via mutation
      parentIndex <- tournament(fitnessValues, popSize, tournamentSize)
      child <- mutationFunction(pop[[parentIndex]])
    }

    childFitness <- fitnessFunction(child)
    bestFitness <- if (childFitness < bestFitness) childFitness else bestFitness

    offspringIndex <- negativeTournament(fitnessValues, popSize, tournamentSize)
    fitnessValues[offspringIndex] <- childFitness
    pop[[offspringIndex]] <- child

    if (archive) {
      archiveList[[length(archiveList) + 1]] <- list(individual = child,
                                                     fitness = childFitness)
    }

    # Apply restart strategy...
    if (restartCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                         evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
      restartResult <- restartStrategy(fitnessFunction, pop, popSize, functionSet, inputVariables, constantSet)
      pop <- restartResult[[1]]
      elite <- joinElites(restartResult[[2]], elite, eliteSize, fitnessFunction)
      logFunction("restart")
    }

    timeElapsed <- proc.time()["elapsed"] - startTime
    stepNumber <- 1 + stepNumber
    evaluationNumber <- 1 + evaluationNumber
    progressMonitor(pop, list(fitnessValues = fitnessValues), fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
  }
  
  elite <- joinElites(pop, elite, eliteSize, fitnessFunction) # insert pop into elite at end of run
  logFunction("Genetic programming evolution run FINISHED after %i evolution steps, %i fitness evaluations and %s.",
              stepNumber, evaluationNumber, formatSeconds(timeElapsed))

  ## Return result list...
  list(timeElapsed = timeElapsed,
       stepNumber = stepNumber,
       evaluationNumber = evaluationNumber,
       bestFitness = bestFitness,
       population = pop,
       fitnessValues = fitnessValues,
       elite = elite,
       archiveList = archiveList,
       searchHeuristicResults = list())
}

##' Comma Evolution Strategy Search Heuristic for RGP
##'
##' The search-heuristic, i.e. the concrete GP search algorithm, is a modular component of RGP.
##' \code{makeCommaEvolutionStrategySearchHeuristic} creates a RGP search-heuristic that implements a
##' (mu, lambda) Evolution Strategy. The lambda parameter is fixed to the population size.
##' TODO description based on Luke09a
##'
##' @param mu The number of surviving parents for the Evolution Strategy search-heuristic. Note that
##'   with \code{makeCommaEvolutionStrategySearchHeuristic}, lambda is fixed to the population size,
##'   i.e. \code{length(pop)}.
##' @return An RGP search heuristic.
##'
##' @export
makeCommaEvolutionStrategySearchHeuristic <- function(mu = 1)
function(logFunction, stopCondition, pop, fitnessFunction,
         mutationFunction, crossoverFunction,
         functionSet, inputVariables, constantSet,
         archive, extinctionPrevention,
         elite, eliteSize,
         restartCondition, restartStrategy,
         breedingFitness, breedingTries,
         progressMonitor) {
  logFunction("STARTING genetic programming evolution run (Evolution Strategy search-heuristic) ...")
  
  ## Tool functions...
  truncationSelect <- function(mu, fitnessValues) order(fitnessValues)[1:mu]

  ## Global variables...
  lambda <- length(pop)
  if (mu > lambda) stop("makeCommaEvolutionStrategySearchHeuristic: mu must be less or equal to the population size")
  childrenPerParent <- ceiling(lambda / mu)
  fitnessValues <- sapply(pop, fitnessFunction)

  ## Initialize statistic counters...
  stepNumber <- 1
  evaluationNumber <- 0
  timeElapsed <- 0
  archiveList <- list() # the archive of all individuals selected in this run, only used if archive == TRUE
  archiveIndexOf <- function(archive, individual)
    Position(function(a) identical(body(a$individual), body(individual)), archive)
  bestFitness <- min(fitnessValues) # best fitness value seen in this run, if multi-criterial, only the first component counts
  startTime <- proc.time()["elapsed"]

  ## Execute GP run...
  while (!stopCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                        evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
    parentIndices <- truncationSelect(mu, fitnessValues)
    nextGeneration <- list()
    for (i in 1:mu) {
      nextGeneration <- c(nextGeneration,
                          replicate(childrenPerParent, mutationFunction(pop[[i]])))
    }
    elite <- joinElites(pop, elite, eliteSize, fitnessFunction) # insert pop into elite
    pop <- nextGeneration[1:lambda] # replace entire population with next generation
    fitnessValues <- sapply(pop, fitnessFunction)
    bestFitness <- if (min(fitnessValues) < bestFitness) min(fitnessValues) else bestFitness
 
    if (archive) {
      stop("not implemented")
    }

    # Apply restart strategy...
    if (restartCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                         evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
      populationSize <- length(pop)
      restartResult <- restartStrategy(fitnessFunction, pop, populationSize, functionSet, inputVariables, constantSet)
      pop <- restartResult[[1]]
      elite <- joinElites(restartResult[[2]], elite, eliteSize, fitnessFunction)
      logFunction("restart")
    }

    timeElapsed <- proc.time()["elapsed"] - startTime
    stepNumber <- 1 + stepNumber
    evaluationNumber <- lambda + evaluationNumber
    progressMonitor(pop, list(fitnessValues = fitnessValues), fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
  }
 
  elite <- joinElites(pop, elite, eliteSize, fitnessFunction) # insert pop into elite at end of run
  logFunction("Genetic programming evolution run FINISHED after %i evolution steps, %i fitness evaluations and %s.",
              stepNumber, evaluationNumber, formatSeconds(timeElapsed))

  ## Return result list...
  list(timeElapsed = timeElapsed,
       stepNumber = stepNumber,
       evaluationNumber = evaluationNumber,
       bestFitness = bestFitness,
       population = pop,
       fitnessValues = fitnessValues,
       elite = elite,
       archiveList = archiveList,
       searchHeuristicResults = list())
}

##' Age Fitness Complexity Pareto GP Search Heuristic for RGP
##'
##' The search-heuristic, i.e. the concrete GP search algorithm, is a modular component of RGP.
##' \code{makeAgeFitnessComplexityParetoGpSearchHeuristic} creates a RGP search-heuristic that implements
##' a generational evolutionary multi objective optimization algorithm (EMOA) that selects on three criteria:
##' Individual age, individual fitness, and individual complexity.
##'
##' @param lambda The number of children to create in each generation (\code{50} by default).
##' @param crossoverProbability The crossover probability for search-heuristics that support
##'   this setting (i.e. TinyGP). Defaults to \code{0.5}.
##' @param enableComplexityCriterion Whether to enable the complexity criterion in multi-criterial
##'   search heuristics.
##' @param enableAgeCriterion Whether to enable the age criterion in multi-criterial search heuristics.
##' @param ndsParentSelectionProbability The probability to use non-dominated sorting to select parents
##'   for each generation. When set to \code{0.0}, parents are selected by uniform random
##'   sampling without replacement every time. Defaults to \code{1.0}.
##' @param ndsSelectionFunction The function to use for non-dominated sorting in Pareto GP selection.
##'   Defaults to \code{nds_cd_selection}.
##' @param complexityMeasure The complexity measure, a function of signature \code{function(ind, fitness)}
##'   returning a single numeric value.
##' @param ageMergeFunction The function used for merging ages of crossover children, defaults
##'   to \code{max}.
##' @param newIndividualsPerGeneration The number of new individuals per generation to
##'   insert into the population. Defaults to \code{50} if \code{enableAgeCriterion == TRUE}
##'   else to \code{0}.
##' @param newIndividualsMaxDepth The maximum depth of new individuals inserted into the
##'   population.
##' @param newIndividualFactory The factory function for creating new individuals. Defaults
##'   to \code{makePopulation}.
##' @return An RGP search heuristic.
##'
##' @export
##' @import emoa
makeAgeFitnessComplexityParetoGpSearchHeuristic <- function(lambda = 50,
                                                            crossoverProbability = 0.5,
                                                            enableComplexityCriterion = TRUE,
                                                            enableAgeCriterion = FALSE,
                                                            ndsParentSelectionProbability = 0.0,
                                                            ndsSelectionFunction = nds_cd_selection,
                                                            complexityMeasure = function(ind, fitness) fastFuncVisitationLength(ind),
                                                            ageMergeFunction = max,
                                                            newIndividualsPerGeneration = if (enableAgeCriterion) 50 else 0,
                                                            newIndividualsMaxDepth = 8,
                                                            newIndividualFactory = makePopulation)
function(logFunction, stopCondition, pop, fitnessFunction,
         mutationFunction, crossoverFunction,
         functionSet, inputVariables, constantSet,
         archive, extinctionPrevention,
         elite, eliteSize,
         restartCondition, restartStrategy,
         breedingFitness, breedingTries,
         progressMonitor) {
  logFunction("STARTING genetic programming evolution run (Age/Fitness/Complexity Pareto GP search-heuristic) ...")

  ## Initialize run-global variables...
  mu <- length(pop)
  if (mu < 2 * lambda) stop("makeAgeFitnessComplexityParetoGpSearchHeuristic: condition mu >= 2 * lambda must be fulfilled")
  fitnessValues <- as.numeric(sapply(pop, fitnessFunction))
  complexityValues <- as.numeric(Map(complexityMeasure, pop, fitnessValues))
  #complexityValues[is.infinite(fitnessValues)] <- Inf # TODO test 
  ageValues <- integer(mu) # initialize ages with zeros
  #ageValues[is.infinite(fitnessValues)] <- Inf # TODO test 
  fastIndividualFactory <- function() {
    makeClosure(.Call("initialize_expression_grow_R",
                      as.list(functionSet$nameStrings),
                      as.integer(functionSet$arities),
                      as.list(inputVariables$nameStrings),
                      -10.0, 10.0,
                      0.8, 0.2,
                      as.integer(newIndividualsMaxDepth)),
                as.list(inputVariables$nameStrings))
  }

  ## Initialize statistic counters...
  stepNumber <- 1
  evaluationNumber <- 0
  timeElapsed <- 0
  archiveList <- list() # the archive of all individuals selected in this run, only used if archive == TRUE
  archiveIndexOf <- function(archive, individual)
    Position(function(a) identical(body(a$individual), body(individual)), archive)
  bestFitness <- min(fitnessValues) # best fitness value seen in this run, if multi-criterial, only the first component counts
  startTime <- proc.time()["elapsed"]

  ## Execute GP run...
  while (!stopCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                        evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {

    # Select 2 * lambda parent individuals...
    parentIndices <- if (runif(1) < ndsParentSelectionProbability) {
      # Select 2 * lambda parent individuals by non-dominated sorting...
      indicesToRemove <- selectIndividualsForReplacement(fitnessValues, complexityValues, ageValues,
                                                         enableComplexityCriterion, enableAgeCriterion,
                                                         ndsSelectionFunction,
                                                         mu - (2 * lambda))
      indicesToKeep <- setdiff(1:mu, indicesToRemove)
      indicesToKeep
    } else {
      # Sample (without replacement) 2 * lambda parent inviduals...
      sample(1:mu, 2 * lambda, replace = FALSE)
    }
    allParentIndices <- 1:(2* lambda)
    motherIndices <- parentIndices[allParentIndices %% 2 == 1] # mothers are odd parent indicies
    fatherIndices <- parentIndices[allParentIndices %% 2 == 0] # fathers are even parent indices

    # Create children individuals...
    children <- Map(function(motherIndex, fatherIndex) {
                      if (runif(1) < crossoverProbability) {
                        # create child via crossover
                        crossoverFunction(pop[[motherIndex]], pop[[fatherIndex]],
                                          breedingFitness = breedingFitness,
                                          breedingTries = breedingTries)
                      } else {
                        # create child via mutation
                        mutationFunction(pop[[motherIndex]])
                      }
                    },
                    motherIndices, fatherIndices)

    # Evaluate children individuals...
    childrenFitnessValues <- as.numeric(sapply(children, fitnessFunction))
    childrenComplexityValues <- as.numeric(Map(complexityMeasure, children, childrenFitnessValues))
    childrenAgeValues <- 1 + as.double(Map(ageMergeFunction, ageValues[motherIndices], ageValues[fatherIndices]))

    # Create and evaluate new individuals...
    newIndividuals <- newIndividualFactory(newIndividualsPerGeneration, functionSet, inputVariables, constantSet,
                                           maxfuncdepth = newIndividualsMaxDepth,
                                           extinctionPrevention = extinctionPrevention,
                                           breedingFitness = breedingFitness, breedingTries = breedingTries,
                                           funcfactory = fastIndividualFactory)
    newIndividualsFitnessValues <- as.numeric(sapply(newIndividuals, fitnessFunction))
    newIndividualsComplexityValues <- as.numeric(Map(complexityMeasure, newIndividuals, newIndividualsFitnessValues))
    newIndividualsAgeValues <- integer(newIndividualsPerGeneration) # initialize ages with zeros

    # Create the pool of individuals to select the next generation from...
    pool <- c(pop, children, newIndividuals)
    poolFitnessValues <- c(fitnessValues, childrenFitnessValues, newIndividualsFitnessValues)
    poolComplexityValues <- c(complexityValues, childrenComplexityValues, newIndividualsComplexityValues) 
    #poolComplexityValues[is.infinite(poolFitnessValues)] <- Inf # TODO test 
    poolAgeValues <- c(ageValues, childrenAgeValues, newIndividualsAgeValues)
    #poolAgeValues[is.infinite(poolFitnessValues)] <- Inf # TODO test 

    # Sort the pool via the non-domination relation and select individuals for removal...
    poolIndicesToRemove <- selectIndividualsForReplacement(poolFitnessValues, poolComplexityValues, poolAgeValues,
                                                           enableComplexityCriterion, enableAgeCriterion,
                                                           ndsSelectionFunction,
                                                           lambda + newIndividualsPerGeneration)

    # Replace current population with next generation...
    pop <- pool[-poolIndicesToRemove]
    fitnessValues <- poolFitnessValues[-poolIndicesToRemove]
    complexityValues <- poolComplexityValues[-poolIndicesToRemove]
    ageValues <- poolAgeValues[-poolIndicesToRemove]
    if (min(fitnessValues) < bestFitness) bestFitness <- min(fitnessValues) # update best fitness

    # Apply restart strategy...
    if (restartCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                         evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
      restartResult <- restartStrategy(fitnessFunction, pop, mu, functionSet, inputVariables, constantSet)
      pop <- restartResult[[1]]
      fitnessValues <- as.numeric(sapply(pop, fitnessFunction))
      complexityValues <- as.numeric(Map(complexityMeasure, pop, fitnessValues))
      ageValues <- integer(mu) # initialize ages with zeros
      logFunction("restart")
    }

    timeElapsed <- proc.time()["elapsed"] - startTime
    stepNumber <- 1 + stepNumber
    evaluationNumber <- lambda + newIndividualsPerGeneration + evaluationNumber
    progressMonitor(pop, list(fitnessValues = fitnessValues, complexityValues = complexityValues, ageValues = ageValues, poolFitnessValues = poolFitnessValues, poolComplexityValues = poolComplexityValues, poolAgeValues = poolAgeValues), fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed, poolIndicesToRemove)
  }
 
  elite <- joinElites(pop, elite, eliteSize, fitnessFunction) # insert pop into elite at end of run
  bestFitness <- min(fitnessValues)
  logFunction("Genetic programming evolution run FINISHED after %i evolution steps, %i fitness evaluations and %s.",
              stepNumber, evaluationNumber, formatSeconds(timeElapsed))

  ## Return result list...
  list(timeElapsed = timeElapsed,
       stepNumber = stepNumber,
       evaluationNumber = evaluationNumber,
       bestFitness = bestFitness,
       population = pop,
       fitnessValues = fitnessValues,
       elite = elite,
       archiveList = archiveList,
       searchHeuristicResults = list(fitnessValues = fitnessValues,
                                     complexityValues = complexityValues,
                                     ageValues = ageValues))
}

##' Archive-based Pareto Tournament Search Heuristic for RGP
##'
##' The search-heuristic, i.e. the concrete GP search algorithm, is a modular component of RGP.
##' \code{makeArchiveBasedParetoTournamentSearchHeuristic} creates a RGP search-heuristic that implements
##' a archive-based Pareto tournament multi objective optimization algorithm (EMOA) that selects on three 
##' criteria: Individual fitness, individual complexity and individual age.
##'
##' @param archiveSize The number of individuals in the archive, defaults to \code{50}.
##' @param popTournamentSize The size of the Pareto tournaments for selecting individuals
##'   for reproduction from the population.
##' @param archiveTournamentSize The size of the Pareto tournaments for selecting individuals
##'   for reproduction from the archive.
##' @param crossoverRate The probabilty to do crossover with an archive member instead of mutation of an
##'   archive member.
##' @param enableComplexityCriterion Whether to enable the complexity criterion in multi-criterial
##'   search heuristics.
##' @param complexityMeasure The complexity measure, a function of signature \code{function(ind, fitness)}
##'   returning a single numeric value.
##' @param ndsSelectionFunction The function to use for non-dominated sorting in Pareto GP selection.
##'   Defaults to \code{nds_cd_selection}.
##' @return An RGP search heuristic.
##'
##' @export
##' @import emoa
makeArchiveBasedParetoTournamentSearchHeuristic <- function(archiveSize = 50,
                                                            popTournamentSize = 5,
                                                            archiveTournamentSize = 3,
                                                            crossoverRate = 0.95,
                                                            enableComplexityCriterion = TRUE,
                                                            complexityMeasure = function(ind, fitness) fastFuncVisitationLength(ind),
                                                            ndsSelectionFunction = nds_cd_selection)
function(logFunction, stopCondition, pop, fitnessFunction,
         mutationFunction, crossoverFunction,
         functionSet, inputVariables, constantSet,
         archive, extinctionPrevention,
         elite, eliteSize,
         restartCondition, restartStrategy,
         breedingFitness, breedingTries,
         progressMonitor) {
  logFunction("STARTING genetic programming evolution run (Archive-based Pareto tournament GP search-heuristic) ...")

  ## Initialize run-global variables...
  mu <- length(pop)
  if (mu < archiveSize) stop("makeAgeFitnessComplexityParetoGpSearchHeuristic: population size (mu) must be larger than or equal to archive size")
  popFitnessValues <- as.numeric(sapply(pop, fitnessFunction))
  popComplexityValues <- as.numeric(Map(complexityMeasure, pop, popFitnessValues))
  popAgeValues <- integer(mu) # initialize ages with zeros

  ## Initialize statistic counters...
  stepNumber <- 1
  evaluationNumber <- 0
  timeElapsed <- 0
  archiveList <- list() # the archive of all individuals selected in this run, only used if archive == TRUE
  archiveIndexOf <- function(archive, individual)
    Position(function(a) identical(body(a$individual), body(individual)), archive)
  bestFitness <- min(popFitnessValues) # best fitness value seen in this run, if multi-criterial, only the first component counts
  startTime <- proc.time()["elapsed"]

  ## Initialize archive with best individuals of population...
  worstPopIndices <- selectIndividualsForReplacement(popFitnessValues, popComplexityValues, popAgeValues,
                                                     enableComplexityCriterion, FALSE, # TODO add age criterion
                                                     ndsSelectionFunction, mu - archiveSize)
  archive <- pop[-worstPopIndices] # initialize archive with best individuals from population
  archiveFitnessValues <- popFitnessValues[-worstPopIndices]
  archiveComplexityValues <- popComplexityValues[-worstPopIndices]
  archiveAgeValues <- popAgeValues[-worstPopIndices]

  allAgeValues <- c(popAgeValues, archiveAgeValues)

  ## Execute GP run...
  while (!stopCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                        evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
    # Generate the next generation's population by crossover and mutation...
    childrenPop <- list()
    for (i in 1:mu) {
      if (runif(1) <= crossoverRate) { # create individual by crossover
        archiveParentIndex <- tournament(archiveFitnessValues, archiveSize, archiveTournamentSize)
        popParentIndex <- tournament(popFitnessValues, mu, popTournamentSize)
        childrenPop[[i]] <- crossoverFunction(archive[[archiveParentIndex]], pop[[popParentIndex]],
                                              breedingFitness = breedingFitness,
                                              breedingTries = breedingTries)
      } else { # create individual by mutation
        parentIndex <- tournament(archiveFitnessValues, archiveSize, archiveTournamentSize)
        childrenPop[[i]] <- mutationFunction(archive[[parentIndex]])
      }
    }

    # Update pop and archive...
    pop <- childrenPop
    popFitnessValues <- as.numeric(sapply(pop, fitnessFunction))
    popComplexityValues <- as.numeric(Map(complexityMeasure, pop, popFitnessValues))
    allIndividuals <- c(pop, archive)
    allFitnessValues <- c(popFitnessValues, archiveFitnessValues)
    allComplexityValues <- c(popComplexityValues, archiveComplexityValues)
    worstIndices <- selectIndividualsForReplacement(allFitnessValues, allComplexityValues, allAgeValues,
                                                    enableComplexityCriterion, FALSE, # TODO add age criterion
                                                    ndsSelectionFunction, mu)
    archive <- allIndividuals[-worstIndices]
    archiveFitnessValues <- allFitnessValues[-worstIndices]
    archiveComplexityValues <- allComplexityValues[-worstIndices]
    if (min(archiveFitnessValues) < bestFitness) bestFitness <- min(archiveFitnessValues) # update best fitness

    # Apply restart strategy...
    if (restartCondition(pop = pop, fitnessFunction = fitnessFunction, stepNumber = stepNumber,
                         evaluationNumber = evaluationNumber, bestFitness = bestFitness, timeElapsed = timeElapsed)) {
      restartResult <- restartStrategy(fitnessFunction, pop, mu, functionSet, inputVariables, constantSet)
      pop <- restartResult[[1]]
      logFunction("restart")
    }

    timeElapsed <- proc.time()["elapsed"] - startTime
    stepNumber <- 1 + stepNumber
    evaluationNumber <- mu + evaluationNumber
    progressMonitor(pop, list(fitnessValues = popFitnessValues, complexityValues = popComplexityValues, ageValues = popAgeValues), fitnessFunction, stepNumber, evaluationNumber, bestFitness, timeElapsed)
  } # evolution main loop
 
  bestFitness <- min(archiveFitnessValues)
  logFunction("Genetic programming evolution run FINISHED after %i evolution steps, %i fitness evaluations and %s.",
              stepNumber, evaluationNumber, formatSeconds(timeElapsed))

  ## Return result list...
  list(timeElapsed = timeElapsed,
       stepNumber = stepNumber,
       evaluationNumber = evaluationNumber,
       bestFitness = bestFitness,
       population = pop,
       fitnessValues = popFitnessValues,
       elite = archive,
       archiveList = archiveList,
       searchHeuristicResults = list(fitnessValues = popFitnessValues,
                                     complexityValues = popComplexityValues))
}


## Tool functions...
randomIndex <- function(maxIndex) as.integer(runif(1, min = 1, max = maxIndex + 1)) 

tournament <- function(fitnessValues, popSize, tournamentSize) {
  bestIndex <- randomIndex(popSize) 
  bestFitness <- Inf
  for (i in 1:tournamentSize) {
    competitorIndex <- randomIndex(popSize)
    if (fitnessValues[competitorIndex] < bestFitness) {
      bestFitness <- fitnessValues[competitorIndex]
      bestIndex <- competitorIndex
    }
  }
  bestIndex
}

negativeTournament <- function(fitnessValues, popSize, tournamentSize) {
  worstIndex <- randomIndex(popSize) 
  worstFitness <- -Inf
  for (i in 1:tournamentSize) {
    competitorIndex <- randomIndex(popSize)
    if (fitnessValues[competitorIndex] > worstFitness) {
      worstFitness <- fitnessValues[competitorIndex]
      worstIndex <- competitorIndex
    }
  }
  worstIndex
}

selectIndividualsForReplacement <- function(fitnessValues, complexityValues, ageValues,
                                            enableComplexityCriterion, enableAgeCriterion,
                                            ndsSelectionFunction,
                                            n) {
  points <- if (enableAgeCriterion & enableComplexityCriterion)
    rbind(fitnessValues, complexityValues, ageValues)
  else if (enableComplexityCriterion)
    rbind(fitnessValues, complexityValues)
  else if (enableAgeCriterion)
    rbind(fitnessValues, ageValues)
  else
    rbind(fitnessValues)
  indicesToRemove <- if (!enableComplexityCriterion & !enableAgeCriterion) {
    # single-criterial case
    order(points, decreasing = TRUE)[1:n]
  } else {
    # multi-criteral case
    ndsSelectionFunction(points, n)
  }
  return (indicesToRemove)
}


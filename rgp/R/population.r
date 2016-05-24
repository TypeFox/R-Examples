## population.R
##   - Functions for handling GP populations
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Classes for populations of individuals represented as functions
##'
##' \code{makePopulation} creates a population of untyped individuals, whereas
##' \code{makeTypedPopulation} creates a population of typed individuals.
##' \code{fastMakePopulation} is a faster variant of \code{makePopulation} with
##'   fewer options.
##' \code{print.population} prints the population.
##' \code{summary.population} returns a summary view of a population.
##'
##' @param size The population size in number of individuals.
##' @param type The (range) type of the individual functions to create.
##' @param funcset The function set.
##' @param inset The set of input variables.
##' @param conset The set of constant factories.
##' @param constMin For \code{fastMakePopulation}, the minimum constant to create.
##' @param constMax For \code{fastMakePopulation}, the maximum constant to create.
##' @param maxfuncdepth The maximum depth of the functions of the new population.
##' @param constprob The probability of generating a constant in a step of growth, if no subtree
##'   is generated. If neither a subtree nor a constant is generated, a randomly chosen input variable
##'   will be generated. Defaults to \code{0.2}.
##' @param breedingFitness A breeding function. See the documentation for
##'   \code{\link{geneticProgramming}} for details.
##' @param breedingTries The number of breeding steps.
##' @param extinctionPrevention When set to \code{TRUE}, initialization will try 
##'   to prevent duplicate individuals from occurring in the population. Defaults to \code{FALSE}, as
##'   this operation might be expensive with larger population sizes.
##' @param funcfactory A factory for creating the functions of the new population.
##'   Defaults to Koza's "ramped half-and-half" initialization strategy.
##' @param x The population to print.
##' @param object The population to summarize.
##' @param ... Additional parameters to the \code{\link{print}} or \code{\link{summary}}
##'   (passed on to their default implementation).
##' @return A new population of functions.
##'
##' @rdname populationCreation
##' @export
makePopulation <- function(size, funcset, inset, conset,
                           maxfuncdepth = 8, constprob = 0.2,
                           breedingFitness = function(individual) TRUE,
                           breedingTries = 50,
                           extinctionPrevention = FALSE,
                           funcfactory = NULL) {
  funcfactory <- if (is.null(funcfactory)) {
    function() randfuncRampedHalfAndHalf(funcset, inset, conset, maxfuncdepth, constprob = constprob,
                                         breedingFitness = breedingFitness, breedingTries = breedingTries)
  } else {
    funcfactory
  }
  pop <- if (size <= 0) {
    list()
  } else if (extinctionPrevention) {
    resultPop <- list()
    l <- 0
    repeat {
      candidate <- funcfactory()
      if (!contains(resultPop, candidate)) {
        resultPop <- c(resultPop, candidate)
        l <- l + 1
      }
      if (l == size) break()
    }
    resultPop
  } else {
    lapply(1:size, function(i) funcfactory())
  }
  class(pop) <- c("untypedPopulation", "population", "list")
  pop
}

##' @rdname populationCreation
##' @export
fastMakePopulation <- function(size, funcset, inset, maxfuncdepth, constMin, constMax) {
  Map(function(i) makeClosure(.Call("initialize_expression_grow_R",
                                    as.list(funcset$nameStrings),
                                    as.integer(funcset$arities),
                                    as.list(inset$nameStrings),
                                    constMin, constMax,
                                    0.8, 0.2,
                                    as.integer(maxfuncdepth)),
                              as.list(inset$nameStrings)), 1:size)
}

##' @rdname populationCreation
##' @export
makeTypedPopulation <- function(size, type, funcset, inset, conset,
                                maxfuncdepth = 8, constprob = 0.2,
                                breedingFitness = function(individual) TRUE,
                                breedingTries = 50,
                                extinctionPrevention = FALSE,
                                funcfactory = NULL) {
  funcfactory <- if (is.null(funcfactory)) {
    function() randfuncTypedRampedHalfAndHalf(type, funcset, inset, conset,
                                              maxfuncdepth, constprob = constprob,
                                              breedingFitness = breedingFitness, breedingTries = breedingTries)
  } else {
    funcfactory
  }
  pop <- makePopulation(size, funcset, inset, conset, maxfuncdepth, funcfactory = funcfactory, constprob = constprob,
                        breedingFitness = breedingFitness, breedingTries = breedingTries,
                        extinctionPrevention = extinctionPrevention)
  class(pop) <- c("typedPopulation", "population", "list")
  pop
}

##' @rdname populationCreation
##' @method print population 
##' @S3method print population 
##' @export
print.population <- function(x, ...) {
  print.default(x, ...) # TODO
}

##' @rdname populationCreation
##' @method summary population 
##' @S3method summary population 
##' @export
summary.population <- function(object, ...) {
  # TODO
  value <- list(class(object)[1], length(object))
  names(value) <- c("Class", "Size")
  class(value) <- "table"
  value
}

##' Fitness/Complexity plot for populations
##'
##' Plots the fitness against the complexity of each individual in a population.
##'
##' @param pop A population to plot.
##' @param fitnessFunction The function to calculate an individual's fitness with.
##' @param complexityFunction The function to calculate an individual's complexity with.
##' @param showIndices Whether to show the population index of each individual.
##' @param showParetoFront Whether to highlight the pareto front in the plot.
##' @param hideOutliers If \code{N = hideOutliers > 0}, hide outliers from the plot using
##'   a "N * IQR" criterion.
##' @param ... Additional parameters for the underlying call to \code{\link{plot}}.
##'
##' @import emoa
##' @export
plotPopulationFitnessComplexity <- function(pop, fitnessFunction,
                                            complexityFunction = fastFuncVisitationLength,
                                            showIndices = TRUE,
                                            showParetoFront = TRUE,
                                            hideOutliers = 0, ...) {
  popFit <- as.vector(lapply(pop, fitnessFunction), mode = "numeric")[1] # only the first fitness component
  popCpx <- as.vector(lapply(pop, complexityFunction), mode = "numeric")
  popFitLim <-
    if (hideOutliers)
      c(0, min(max(popFit), quantile(popFit, 0.75) + hideOutliers * IQR(popFit)))
    else NULL
  popCpxLim <-
    if (hideOutliers)
      c(0, min(max(popCpx), quantile(popCpx, 0.75) + hideOutliers * IQR(popCpx)))
    else NULL
  popNds <- nds_rank(rbind(popFit, popCpx))
  points <- cbind(popFit, popCpx)
  colnames(points) <- c("Fitness", "Complexity")
  popPch <- rep(1, length(pop))
  if (showParetoFront) stop("showParetoFront not implemented") # TODO
  if (showParetoFront) popPch <- replace(popPch, which(popNds == 1), 16)
  plot(points, pch = popPch, xlim = popFitLim, ylim = popCpxLim, ...)
  if (showIndices) text(x = popFit, y = popCpx, pos = 1, cex = 0.6,
                        xlim = popFitLim, ylim = popCpxLim, ...)
}

##' Calculate the fitness value of each individual in a population
##'
##' @param pop A population of functions.
##' @param fitnessfunc The fitness function.
##' @return A list of fitness function values in the same order as \code{pop}.
##'
##' @export
popfitness <- function(pop, fitnessfunc) sapply(pop, fitnessfunc)

## breeding.R
##   - Functions for "breeding" individuals by repeated application of init-
##     ialization or variation operators
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Breeding of GP individuals
##'
##' Breeds GP individuals by repeated application of an individual factory function.
##' \code{individualFactory}. The \code{breedingFitness} must be a function of domain
##' logical (a single boolean value) or numeric (a single real number). In case of
##' a boolean breeding function, candidate individuals are created via the
##' \code{individualFactory} function and tested by the \code{breedingFitness} predicate
##' until the \code{breedingFitness} predicate is \code{TRUE} or \code{breedingTries} tries
##' were done, in which case the last individual created and tested is returned. In case
##' of a numerical breeding function, \code{breedingTries} individuals are created and
##' evaluated by the \code{breedingFitness} function. The individual with the minimal
##' breeding fitness is returned.
##'
##' @param individualFactory A function of no parameters that returns a single GP
##'   individual.
##' @param breedingFitness Either a function that takes a GP individual as its only
##'   parameter and returns a single logical value or a function that takes a GP
##'   individual as its only parameter and returns a single real value.
##' @param breedingTries The number of breeding steps to perform. In case of a
##'   boolean \code{breedingFitness} function, the actual number of breeding
##'   steps performed may be lower then this number (see the details).
##' @param warnOnFailure Whether to issue a warning when a boolean \code{breedingFitness}
##'   predicate was not fulfilled after \code{breedingTries} tries.
##' @param stopOnFailure Whether to stop with an error message when a boolean
##'   \code{breedingFitness} predicate was not fulfilled after \code{breedingTries} tries.
##' @return The GP individual that was bred.
##'
##' @rdname breeding
##' @export
breed <- function(individualFactory, breedingFitness, breedingTries,
                  warnOnFailure = TRUE, stopOnFailure = FALSE) {
  candidateIndividual <- individualFactory()
  candidateFitness <- breedingFitness(candidateIndividual)
  selectedIndividual <- candidateIndividual
  repeats <- 1
  
  if (is.logical(candidateFitness)) {
    while (!candidateFitness && repeats < breedingTries) {
      candidateIndividual <- individualFactory()
      candidateFitness <- breedingFitness(candidateIndividual)
      selectedIndividual <- candidateIndividual
      repeats <- repeats + 1
    }
    if (!candidateFitness && (warnOnFailure || stopOnFailure)) {
      if (!stopOnFailure)
        warning("breed: Boolean breeding fitness predicate unfulfilled after ", breedingTries, " tries.")
      else
        stop("breed: Boolean breeding fitness predicate unfulfilled after ", breedingTries, " tries.")
    }
  } else if (is.numeric(candidateFitness)) {
    bestFitness <- candidateFitness
    while (repeats < breedingTries) {
      candidateIndividual <- individualFactory()
      candidateFitness <- breedingFitness(candidateIndividual)
      if (candidateFitness < bestFitness) {
        selectedIndividual <- candidateIndividual
        bestFitness <- candidateFitness
      }
      repeats <- repeats + 1
    }
  } else {
    stop("breed: BreedingFitness must be a function that returns a logical or a numeric.")
  }
  
  selectedIndividual
}

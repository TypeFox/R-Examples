## selection.R
##   - GP selection strategies
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' GP selection functions
##'
##' A GP selection function determines which individuals in a population should
##' survive, i.e. are selected for variation or cloning, and which individuals
##' of a population should be replaced. Single-objective selection functions base
##' their selection decision on scalar fitness function, whereas multi-objective
##' selection functions support vector-valued fitness functions.
##' Every selection function takes a population and a (possibly vector-valued) fitness
##' function as required arguments. It returns a list of two tables \code{selected}
##' and \code{discarded}, with columns \code{index} and \code{fitness} each. The
##' returned list also contains a single integer \code{numberOfFitnessEvaluations}
##' that contains the number of fitness evaluations used to make the selection (Note
##' that in the multi-objective case, evaluating all fitness functions once counts
##' as a single evaluation). The first table contains the population indices of the
##' individuals selected as survivors, the second table contains the population indices
##' of the individuals that should be discarded and replaced. This definition simplifies
##' the implementation of \emph{steady-state} evolutionary strategies where most of the
##' individuals in a population are unchanged in each selection step. In a GP context,
##' steady-state strategies are often more efficient than generational strategies. 
##'
##' \code{makeTournamentSelection} returns a classic single-objective tournament selection
##'   function.
##' \code{makeMultiObjectiveTournamentSelection} returns a multi-objective tournament selection
##'   function that selects individuals based on multiple objectives.
##' \code{makeComplexityTournamentSelection} returns a multi-objective selection function that
##'   implements the common case of dual-objective tournament selection with high solution
##'   quality as the first objective and low solution complexity as the second objective.
##'
##' @param complexityMeasure The function used to measure the complexity of an individual.
##' @param tournamentSize The number of individuals to randomly select to form a
##'   tournament, defaults to 10 in the single-objective case, 30 in the multi-objective case.
##' @param selectionSize The number of individuals to return as selected.
##' @param tournamentDeterminism The propability \emph{p} for selecting the best individual
##'   in a tournament, must be in the interval (0.0, 1.0]. The best individual is selected
##'   with propability \emph{p}, the second best individual is selected with propability
##'   \emph{p * (1 - p)}, the third best individual ist selected with propability
##'   \emph{p * (1 - p)^2}, and so on. Note that setting \code{tournamentDeterminism}
##'   to \code{1.0} (the default) yields determistic behavior.
##' @param vectorizedFitness If \code{TRUE}, the fitness function is expected to take
##'   a list of individuals as input and return a list of (possible vector-valued) fitnesses
##'   as output.
##' @param rankingStrategy The strategy used to rank individuals based on multiple objectives.
##'   This function must turn a fitness vector (one point per column) into an ordering
##'   permutation (similar to the one returned by \code{order}). Defaults to
##'   \code{orderByParetoCrowdingDistance}.
##' @return A selection function.
##'
##' @rdname selectionFunctions
##' @export
makeTournamentSelection <- function(tournamentSize = 10,
                                    selectionSize = ceiling(tournamentSize / 2),
                                    tournamentDeterminism = 1.0,
                                    vectorizedFitness = FALSE)
  function(population, fitnessFunction) {
    poolIdxs <- sample(length(population), tournamentSize)
    poolFitness <-
      if (vectorizedFitness) {
        fitnessFunction(population[poolIdxs])
      } else {
        sapply(population[poolIdxs], fitnessFunction)
      }
    poolFitness <- sapply(poolFitness, function(fitness) fitness[1]) # only use the first fitness component
    idxFitTable <- cbind(poolIdxs, poolFitness)
    colnames(idxFitTable) <- c("index", "fitness")
    # Sort by (single-objective) fitness...
    sortedIdxFitTable <- idxFitTable[order(idxFitTable[,"fitness"]),]
    # ...then shuffle the ranking depending on tournamentDeterminism:
    shuffledSortedIdxFitTable <-
      sortedIdxFitTable[inversePermutation(nondeterministicRanking(tournamentSize)),]
    # The first selectionSize individuals are selected, the rest are discarded:
    list(selected = matrix(shuffledSortedIdxFitTable[1:selectionSize,], nrow = selectionSize),
         discarded = matrix(shuffledSortedIdxFitTable[-(1:selectionSize),], nrow = tournamentSize - selectionSize),
         numberOfFitnessEvaluations = tournamentSize)
  }

##' @rdname selectionFunctions
##' @export
makeMultiObjectiveTournamentSelection <- function(tournamentSize = 30,
                                                  selectionSize = ceiling(tournamentSize / 2),
                                                  tournamentDeterminism = 1.0,
                                                  vectorizedFitness = FALSE,
                                                  rankingStrategy = orderByParetoCrowdingDistance)
  function(population, fitnessFunction) {
    poolIdxs <- sample(length(population), tournamentSize)
    poolFitness <-
      if (vectorizedFitness) {
        fitnessFunction(population[poolIdxs])
      } else {
        sapply(population[poolIdxs], fitnessFunction) # all fitness components
      }
    if (!is.matrix(poolFitness)) poolFitness <- matrix(poolFitness, ncol = length(poolFitness))
    idxFitTable <- cbind(poolIdxs, t(poolFitness))
    # Sort by (multi-objective) fitness...
    sortedIdxFitTable <- idxFitTable[rankingStrategy(poolFitness),]
    # ...then shuffle the ranking depending on tournamentDeterminism:
    shuffledSortedIdxFitTable <-
      sortedIdxFitTable[inversePermutation(nondeterministicRanking(tournamentSize)),]
    # The first selectionSize individuals are selected, the rest are discarded:
    list(selected = matrix(shuffledSortedIdxFitTable[1:selectionSize,], nrow = selectionSize),
         discarded = matrix(shuffledSortedIdxFitTable[-(1:selectionSize),], nrow = tournamentSize - selectionSize),
         numberOfFitnessEvaluations = tournamentSize)
  }

##' @rdname selectionFunctions
##' @export
makeComplexityTournamentSelection <- function(tournamentSize = 30,
                                              selectionSize = ceiling(tournamentSize / 2),
                                              tournamentDeterminism = 1.0,
                                              vectorizedFitness = FALSE,
                                              rankingStrategy = orderByParetoCrowdingDistance,
                                              complexityMeasure = fastFuncVisitationLength) {
  selectionFunction <- makeMultiObjectiveTournamentSelection(tournamentSize, selectionSize,
                                                             tournamentDeterminism, rankingStrategy)
  # TODO add support for vectorizedFitness
  if (vectorizedFitness) stop("makeComplexityTournamentSelection: vectorizedFitness not implemented")
  function(population, fitnessFunction) {
    complexityAugmentedFintessFunction <- function(individual)
      c(fitnessFunction(individual), complexityMeasure(individual))
    selectionFunction(population, complexityAugmentedFintessFunction)
  }
}

##' Rearrange points via Pareto-based rankings
##'
##' Returns a permutation that rearranges points, given as columns in a value matrix, via
##' Pareto-based ranking. Points are ranked by their Pareto front number. In
##' \code{orderByParetoCrowdingDistance}, ties are then broken by crowding distance,
##' in \code{orderByParetoHypervolumeContribution}, ties are broken by hypervolume
##' contribution.
##'
##' @param values The value matrix to return the ordering permutation for. Each column
##'   represents a point, each row a dimension.
##' @return A permutation to rearrange \code{values} based on a Pareto based ranking.
##'
##' @rdname paretoSorting
##' @import emoa
##' @export
orderByParetoCrowdingDistance <- function(values)
  orderByParetoMeasure(values, measure = crowding_distance)

##' @rdname paretoSorting
##' @import emoa
##' @export
orderByParetoHypervolumeContribution <- function(values)
  orderByParetoMeasure(values, measure = hypervolume_contribution)

##' Rearrange points via an arbitrary Pareto-based ranking
##'
##' Returns a permutation that rearranges points, given as columns in a value matrix, via
##' Pareto-based ranking. Points are ranked by their Pareto front number, ties are broken
##' by the values of \code{measure}.
##'
##' @param values The value matrix to return the ordering permutation for. Each column
##'   represents a point, each row a dimension.
##' @param measure The measure used for ranking points that lie on the same Pareto front,
##'   defaults to \code{crowding_distance}.
##' @return A permutation to rearrange \code{values} based on a Pareto based ranking.
##'
##' @import emoa
orderByParetoMeasure <- function(values, measure = crowding_distance) {
  ndsRanks <- nds_rank(values)
  indicesOrderedByMeasure <- c()
  for (ndsRank in 1:max(ndsRanks)) {
    indicesOfNdsRank <- which(ndsRanks == ndsRank)
    valuesOfNdsRank <- values[,ndsRanks == ndsRank]
    measureRanks <- rank(measure(as.matrix(valuesOfNdsRank)), ties.method = "random")
    indicesOfNdsRankOrderedByMeasure <- indicesOfNdsRank[inversePermutation(measureRanks)]
    indicesOrderedByMeasure <- append(indicesOrderedByMeasure, indicesOfNdsRankOrderedByMeasure)
  }
  indicesOrderedByMeasure
}

##' Create a nondeterministic ranking
##'
##' Create a permutation of the sequence \code{s} = \code{1:l} representing a ranking.
##' If \code{p} = 1, the ranking will be completely deterministic, i.e. equal to
##' \code{1:l}. If \code{p} = 0, the ranking will be completely random. If
##' 0 < \code{p} < 1, the places in the ranking will be determined by iterative
##' weighted sampling without replacement from the sequence \code{s} := \code{1:l}.
##' At each step of this iterated weighted sampling, the first remaining element of
##' \code{s} will be selected with probability \code{p}, the second element with
##' probability \code{p * (1 - p)}, the third element with probability
##' \code{p * (1 - p) ^ 2}, and so forth.
##'
##' @param l The numer of elements in the ranking.
##' @param p The "degree of determinism" of the ranking to create.
##' @return A ranking permutation of the values \code{1:l}.
nondeterministicRanking <- function(l, p = 1) {
  if (p == 0)
    sample(1:l, replace = FALSE) # completely random ranking
  else if (p == 1)
    1:l # completely deterministic ranking
  else {
    # The following more straight-forward code suffers from numeric instability for l >~ 50:
    # pmf <- function(rdet, rl)
    #   if (rl <= 1) rdet else c(rdet, sapply(1:(rl - 1), function(i) rdet * (1 - rdet) ^ i))
    # sample(1:l, replace = FALSE, prob = pmf(p, l)) # TODO numeric bug for l >> 10!
    # This is why we do it the stupid way instead:
    indexPmf <- function(rdet, l) {
      randomPool <- runif(l); randomPool[l] <- -Inf # at least the last index is selected
      which(randomPool <= rdet)[1]
    }
    urn <- 1:l
    ranking <- c()
    for (i in 1:l) {
      sampledIndex <- indexPmf(p, l - i + 1) # Pull an index from the urn...
      ranking <- c(ranking, urn[sampledIndex]) # ...and prepend it to the ranking,...
      urn <- urn[-sampledIndex] # ...then remove it from the urn.
    }
    ranking
  }
}

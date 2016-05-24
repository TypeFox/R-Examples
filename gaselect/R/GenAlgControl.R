#' Control class for the genetic algorithm
#'
#' This class controls the general setup of the genetic algorithm
#' @slot populationSize The number of "chromosomes" in the population (between 1 and 2^16).
#' @slot numGenerations The number of generations to produce (between 1 and 2^16).
#' @slot minVariables The minimum number of variables in the variable subset (between 0 and p - 1 where p is the total number of variables).
#' @slot maxVariables The maximum number of variables in the variable subset (between 1 and p, and greater than \code{minVariables}).
#' @slot elitism The number of absolute best chromosomes to keep across all generations (between 1 and min(\code{populationSize} * \code{numGenerations}, 2^16)).
#' @slot mutationProbability The probability of mutation (between 0 and 1).
#' @slot badSolutionThreshold The child must not be more than \code{badSolutionThreshold} percent worse than the worse parent. If less than 0, the child must be even better than the worst parent.
#' @slot crossover The crossover method to use
#' @slot crossoverId The numeric ID of the crossover method to use
#' @slot maxDuplicateEliminationTries The maximum number of tries to eliminate duplicates
#' @slot verbosity The level of verbosity. 0 means no output at all, 2 is very verbose.
#' @aliases GenAlgControl
#' @rdname GenAlgControl-class
setClass("GenAlgControl", representation(
	populationSize = "integer",
	numGenerations = "integer",
	minVariables = "integer",
	maxVariables = "integer",
	elitism = "integer",
	mutationProbability = "numeric",
	crossover = "character",
	crossoverId = "integer",
	fitnessScaling = "character",
	fitnessScalingId = "integer",
	badSolutionThreshold = "numeric",
	maxDuplicateEliminationTries = "integer",
	verbosity = "integer"
), validity = function(object) {
	errors <- character(0);
	MAXUINT16 <- 2^16; # unsigned 16bit integers are used (uint16_t) in the C++ code

	## Type checks:
	if(object@populationSize < 0L || object@populationSize > MAXUINT16) {
		errors <- c(errors, paste("The population size must be between 0 and", MAXUINT16));
	}

	if(object@numGenerations < 0L || object@numGenerations > MAXUINT16) {
		errors <- c(errors, paste("The number of generations must be between 0 and", MAXUINT16));
	}

	if(object@elitism < 0L || object@elitism > MAXUINT16) {
		errors <- c(errors, paste("'elitism' must be between 0 and", MAXUINT16));
	}

	if(object@minVariables < 0L || object@minVariables > MAXUINT16) {
		errors <- c(errors, paste("The minimal number of variables must be between 0 and", MAXUINT16));
	}

	if(object@maxVariables < 0L || object@maxVariables > MAXUINT16) {
		errors <- c(errors, paste("The maximum number of variables must be strictly larger than the minimum number of variables and between 0 and", MAXUINT16));
	}

	## Sanity checks:

	if(object@populationSize < object@elitism) {
		errors <- c(errors, "The population size must be at least as large as the number of elite solutions");
	}

	if(object@minVariables >= object@maxVariables) {
		errors <- c(errors, "The minimal number of variables must be strictly less than the maximum number");
	}

	if(object@elitism >= object@populationSize * object@numGenerations) {
		errors <- c(errors, "Requested more elite solutions than possible");
	}

	if(object@mutationProbability < 0 || object@mutationProbability >= 1) {
		errors <- c(errors, "The mutation probability must be between 0 and 1 (excluded).");
	}

	if(object@maxDuplicateEliminationTries < 0L) {
		errors <- c(errors, "The maximum number of tries to eliminate duplicates must be greater or equal 0");
	}

	if(object@verbosity < 0L || object@verbosity > 5L) {
		errors <- c(errors, "The verbosity level can not be less than 0 or greater than 5");
	}

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
});


#' Set control arguments for the genetic algorithm
#'
#' The population must be large enough to allow the algorithm to explore the whole solution space. If
#' the initial population is not diverse enough, the chance to find the global optimum is very small.
#' Thus the more variables to choose from, the larger the population has to be.
#'
#' The initial population is generated randomly. Every chromosome uses between \code{minVariables} and
#' \code{maxVariables} (uniformly distributed).
#'
#' If the mutation probability (\code{mutationProbability} is greater than 0, a random number of
#' variables is added/removed according to a truncated geometric distribution to each offspring-chromosome.
#' The resulting distribution of the total number of variables in the subset is not uniform anymore, but almost (the smaller the
#' mutation probability, the more "uniform" the distribution). This should not be a problem for most
#' applications.
#'
#' The user can choose between \code{single} and \code{random} crossover for the mating process. If single crossover
#' is used, a single position is randomly chosen that marks the position to split both parent chromosomes. The child
#' chromosomes are than the concatenated chromosomes from the 1st part of the 1st parent and the 2nd part of the
#' 2nd parent resp. the 2nd part of the 1st parent and the 1st part of the 2nd parent.
#' Random crossover is that a random number of random positions are drawn and these positions are transferred
#' from one parent to the other in order to generate the children.
#'
#' Elitism is a method of enhancing the GA by keeping track of very good solutions. The parameter \code{elitism}
#' specifies how many "very good" solutions should be kept.
#'
#' Before the selection probabilities are determined, the fitness values \eqn{f} of the chromosomes are
#' standardized to the z-scores (\eqn{z = (f - mu) / sd}). Scaling the fitness values afterwards with
#' the exponential function can help the algorithm to faster find good solutions. When setting
#' \code{fitnessScaling} to \code{"exp"}, the (standardized) fitness \eqn{z} will be scaled by \eqn{exp(z)}.
#' This promotes good solutions to get an even higher selection probability, while bad solutions
#' will get an even lower selection probability.
#'
#' @param populationSize The number of "chromosomes" in the population (between 1 and 2^16)
#' @param numGenerations The number of generations to produce (between 1 and 2^16)
#' @param minVariables The minimum number of variables in the variable subset (between 0 and p - 1 where p is the total number of variables)
#' @param maxVariables The maximum number of variables in the variable subset (between 1 and p, and greater than \code{minVariables})
#' @param elitism The number of absolute best chromosomes to keep across all generations (between 1 and min(\code{populationSize} * \code{numGenerations}, 2^16))
#' @param mutationProbability The probability of mutation (between 0 and 1)
#' @param crossover The crossover type to use during mating (see details). Partial matching is performed
#' @param badSolutionThreshold The worst child must not be more than \code{badSolutionThreshold} times worse than the worse parent.
#'			If less than 0, the child must be even better than the worst parent. If the algorithm can't find a better child
#'			in a long time it issues a warning and uses the last found child to continue.
#' @param maxDuplicateEliminationTries The maximum number of tries to eliminate duplicates
#'        (a value of \code{0} or \code{NULL} means that no checks for duplicates are done.
#' @param verbosity The level of verbosity. 0 means no output at all, 2 is very verbose.
#' @param fitnessScaling How the fitness values are internally scaled before the selection probabilities are assigned
#'          to the chromosomes. See the details for possible values and their meaning.
#' @return An object of type \code{\link{GenAlgControl}}
#' @export
#' @example examples/genAlg.R
#' @rdname GenAlgControl-constructor
genAlgControl <- function(populationSize, numGenerations, minVariables, maxVariables,
							elitism = 10L, mutationProbability = 0.01, crossover = c("single", "random"),
							maxDuplicateEliminationTries = 0L, verbosity = 0L, badSolutionThreshold = 2,
							fitnessScaling = c("none", "exp")) {
	if(is.numeric(populationSize)) {
		populationSize <- as.integer(populationSize);
	}

	if(is.numeric(numGenerations)) {
		numGenerations <- as.integer(numGenerations);
	}

	if(is.numeric(minVariables)) {
		minVariables <- as.integer(minVariables);
	}

	if(is.numeric(maxVariables)) {
		maxVariables <- as.integer(maxVariables);
	}

	if(is.numeric(elitism)) {
		elitism <- as.integer(elitism);
	}

	if(is.numeric(verbosity)) {
		verbosity <- as.integer(verbosity);
	}

    if(is.null(maxDuplicateEliminationTries)) {
        maxDuplicateEliminationTries <- 0L;
    }

	crossover <- match.arg(crossover);

	crossoverId <- switch(crossover,
		single = 0L,
		random = 1L
	);

    fitnessScaling <- match.arg(fitnessScaling);
    fitnessScalingId <- switch(fitnessScaling,
		none = 0L,
		exp = 1L
	);

	return(new("GenAlgControl",
				populationSize = populationSize,
				numGenerations = numGenerations,
				minVariables = minVariables,
				maxVariables = maxVariables,
				elitism = elitism,
				mutationProbability = mutationProbability,
				crossover = crossover,
				crossoverId = crossoverId,
				maxDuplicateEliminationTries = as.integer(maxDuplicateEliminationTries),
				badSolutionThreshold = badSolutionThreshold,
				fitnessScaling = fitnessScaling,
				fitnessScalingId = fitnessScalingId,
				verbosity = verbosity));
};

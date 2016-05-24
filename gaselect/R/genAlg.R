#' Result of a genetic algorithm run
#'
#' Return object of a run of the genetic algorithm genAlg
#' @slot subsets Logical matrix with one variable subset per column. The columns are ordered according to their fitness (first column contains the fittest variable-subset).
#' @slot rawFitness Numeric vector with the raw fitness of the corresponding variable subset returned by the evaluator.
#' @slot response The original response vector.
#' @slot covariates The original covariates matrix.
#' @slot evaluator The evaluator used in the genetic algorithm.
#' @slot control The control object.
#' @slot segmentation The segments used by the evaluator. Empty list if the evaluator doesn't use segmentation.
#' @slot seed The seed the algorithm is started with.
#' @aliases GenAlg
#' @include Evaluator.R GenAlgControl.R
#' @import methods
#' @rdname GenAlg-class
setClass("GenAlg", representation(
	subsets = "matrix",
	rawFitness = "numeric",
	rawFitnessEvolution = "matrix",
	response = "numeric",
	covariates = "matrix",
	evaluator = "GenAlgEvaluator",
	control = "GenAlgControl",
	segmentation = "list",
	seed = "integer"
), prototype(
	subsets = matrix(),
	rawFitness = NA_real_
), validity = function(object) {
	errors <- character(0);
	if(!is.numeric(object@response) || !(is.vector(object@response) || is.matrix(object@response) && ncol(object@response) == 1)) {
		errors <- c(errors, "The response variable must be a numeric vector");
	}

	if(!is.numeric(object@covariates) || !is.matrix(object@covariates)) {
		errors <- c(errors, "The covariates must be a numerical matrix");
	}

	if(length(object@response) != nrow(object@covariates)) {
		errors <- c(errors, "The response and the covariates must have the same number of observations");
	}

	if(object@control@minVariables >= ncol(object@covariates)) {
		errors <- c(errors, "The minimum number of variables must be strictly less than the number of available variables");
	}

	if(object@control@maxVariables > ncol(object@covariates)) {
		errors <- c(errors, "The maximum number of variables must be less or equal than the number of available variables");
	}

	dataErrors <- validData(object@evaluator, object)
	if(!is.logical(dataErrors)) {
		errors <- c(errors, dataErrors);
	} else if(dataErrors == FALSE) {
		errors <- c(errors, "The evaluator can not handle this kind of data");
	}

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
});

#' Genetic algorithm for variable subset selection
#'
#' A genetic algorithm to find "good" variable subsets based on internal PLS evaluation or a user specified
#' evaluation function
#'
#' The GA generates an initial "population" of \code{populationSize} chromosomes where each initial
#' chromosome has a random number of randomly selected variables. The fitness of every chromosome is evaluated by
#' the specified evaluator. The default built-in PLS evaluator (see \code{\link{evaluatorPLS}}) is the preferred
#' evaluator.
#' Chromosomes with higher fitness have higher probability of mating with another chromosome. \code{populationSize / 2} couples each create
#' 2 children. The children are created by randomly mixing the parents' variables. These children make up the new generation and are again
#' selected for mating based on their fitness. A total of \code{numGenerations} generations are built this way.
#' The algorithm returns the last generation as well as the best \code{elitism} chromosomes from all generations.
#'
#' @param y The numeric response vector of length n
#' @param X A n x p numeric matrix with all p covariates
#' @param control Options for controlling the genetic algorithm. See \code{\link{genAlgControl}} for details.
#' @param evaluator The evaluator used to evaluate the fitness of a variable subset. See
#'      \code{\link{evaluatorPLS}}, \code{\link{evaluatorLM}} or \code{\link{evaluatorUserFunction}} for details.
#' @param seed Integer with the seed for the random number generator or NULL to automatically seed the RNG
#' @export
#' @import Rcpp
#' @include Evaluator.R GenAlgControl.R formatSegmentation.R
#' @return An object of type \code{\link{GenAlg}}
#' @rdname GenAlg-constructor
#' @useDynLib gaselect
#' @example examples/genAlg.R
genAlg <- function(y, X, control, evaluator = evaluatorPLS(), seed) {
    seed <- as.integer(seed)[1];

    if (!is.numeric(seed) | is.na(seed)) {
        stop("`seed` must be an integer.");
    }

	ret <- new("GenAlg",
		response = y,
		covariates = X,
		evaluator = evaluator,
		control = control,
		seed = seed
	);

	possSubsetCutoff <- 0.85;
	numPossibleSubsets <- sum(choose(ncol(ret@covariates), seq.int(ret@control@minVariables, ret@control@maxVariables)));

    seed <- as.integer(seed);

	if(ret@control@populationSize > possSubsetCutoff * numPossibleSubsets) {
		stop(paste("Requested a population that is almost as large as the number of all possible subsets. The population size can be at most ",
			floor(possSubsetCutoff * 100),
			" of the total number of possible subsets (i.e., ",
			floor(possSubsetCutoff * numPossibleSubsets), ").", sep = ""));
	}

	ctrlArg <- c(toCControlList(ret@control), toCControlList(ret@evaluator));
	ctrlArg$chromosomeSize = ncol(ret@covariates);

	ctrlArg$userEvalFunction <- getEvalFun(ret@evaluator, ret);

	if(ctrlArg$evaluatorClass == 0) {
		res <- .Call("genAlgPLS", ctrlArg, NULL, NULL, seed, PACKAGE = "gaselect");
	} else {
		res <- .Call("genAlgPLS", ctrlArg, ret@covariates, as.matrix(ret@response), ret@seed, PACKAGE = "gaselect");
	}

	ret@subsets <- res$subsets;
	ret@segmentation <- formatSegmentation(ret@evaluator, res$segmentation);
	ret@rawFitness <- res$fitness;
	ret@rawFitnessEvolution <- matrix(res$fitnessEvolution, ncol = 3L, byrow = TRUE, dimnames = list(NULL, c("best", "mean", "std.dev")));

	return(ret);
}

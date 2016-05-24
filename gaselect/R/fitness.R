#' Get the evolution of the fitness
#'
#' Get the fitness of the best / average chromosomes after each generation
#'
#' Returns the progress of the fitness of the best or average chromosome.
#'
#' @param object The \code{\link{GenAlg}} object returned by \code{\link{genAlg}}
#' @param what can be one ore more of \code{"best"} (to return the fitness of the best chromosome for each generation),
#' \code{"mean"} (to return the arithmetic mean fitness during each generation), and \code{"std.dev"} (for
#' the standard deviation of the fitness values in each generation).
#' @param type one of \code{"true"} or \code{"raw"}. \emph{raw} means the raw fitness value used
#' within the GA, while \emph{true} tries to convert it to the standard error of prediction (like
#' \code{\link{fitness}}). If the standard deviation (\code{what = "std.dev"}) is requested, the
#' \code{type} will always be \emph{raw}.
#' @return A vector with the best or average fitness value after each generation
#' @example examples/fitness.R
#' @export
fitnessEvolution <- function(object, what = c("mean", "best", "std.dev"), type = c("true", "raw")) {
    what <- match.arg(what, c("mean", "best", "std.dev"), several.ok = TRUE);
    type <- match.arg(type);
    fit <- object@rawFitnessEvolution[ , what, drop = FALSE];

    if ("std.dev" %in% what || type == "raw") {
        return(fit);
    } else {
        return(apply(fit, 2, function(x) { trueFitnessVal(object@evaluator, x); }));
    }
}

#' Get the fitness of a variable subset
#'
#' Get the internal fitness for all variable subsets
#'
#' This method is used to get the fitness of all variable subsets
#' found by the genetic algorithm.
#'
#' @param object The \code{\link{GenAlg}} object returned by \code{\link{genAlg}}
#' @return A vector with the estimated fitness for each solution
#' @export
#' @example examples/fitness.R
fitness <- function(object) {
    return(trueFitnessVal(object@evaluator, object@rawFitness));
}

#' Get the transformed fitness values
#'
#' Transform the given fitness values according tho the GenAlgEvaluator class
#'
#' This method is used to calculate the true fitness given the GenAlgEvaluator class (as they use
#' different internal fitness measures)
#'
#' @param object The used evaluator, an object with type or with a subtype of \code{\link{GenAlgEvaluator}}
#' @param fitness A numeric vector of fitnesses
#' @return A vector with the true fitness values
#' @docType methods
#' @rdname trueFitnessVal-methods
setGeneric("trueFitnessVal", function(object, fitness) { standardGeneric("trueFitnessVal"); });

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgPLSEvaluator", fitness = "numeric"), function(object, fitness) {
	return(-fitness);
});

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgUserEvaluator", fitness = "numeric"), function(object, fitness) {
	return(object@sepFunction(fitness));
});

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgLMEvaluator", fitness = "numeric"), function(object, fitness) {
	return(fitness);
});

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgFitEvaluator", fitness = "numeric"), function(object, fitness) {
	return(fitness);
});


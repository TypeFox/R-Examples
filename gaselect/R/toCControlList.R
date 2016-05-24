#' Transform the object to a list
#'
#' Get the control list for the C++ procedure genAlgPLS from the object
#'
#' @param object The object
#' @return A list with all items expected by the C++ code
#' @docType methods
#' @include Evaluator.R GenAlgControl.R
#' @rdname toCControlList-methods
setGeneric("toCControlList", function(object) { standardGeneric("toCControlList"); });

#' @rdname toCControlList-methods
setMethod("toCControlList", signature(object = "GenAlgPLSEvaluator"), function(object) {
	return(list(
		"evaluatorClass" = 1L,
		"numReplications" = object@numReplications,
		"innerSegments" = object@innerSegments,
        "outerSegments" = object@outerSegments,
	    "testSetSize" = object@testSetSize,
	    "sdfact" = object@sdfact,
		"plsMethod" = object@methodId,
		"numThreads" = object@numThreads,
        "maxNComp" = object@maxNComp,
		"userEvalFunction" = function() {NULL;},
		"statistic" = 0L
	));
});
#' @rdname toCControlList-methods
setMethod("toCControlList", signature(object = "GenAlgFitEvaluator"), function(object) {
	return(list(
		"evaluatorClass" = 3L,
		"numReplications" = 0L,
		"innerSegments" = object@numSegments,
        "outerSegments" = 0L,
	    "testSetSize" = 0.0,
	    "sdfact" = object@sdfact,
		"plsMethod" = 0L,
		"numThreads" = object@numThreads,
        "maxNComp" = object@maxNComp,
		"userEvalFunction" = function() {NULL;},
		"statistic" = object@statisticId
	));
});

#' @rdname toCControlList-methods
setMethod("toCControlList", signature(object = "GenAlgUserEvaluator"), function(object) {
	return(list(
		"evaluatorClass" = 0L,
		"numReplications" = 0L,
	    "innerSegments" = 0L,
	    "outerSegments" = 0L,
	    "testSetSize" = 0.0,
	    "sdfact" = 0.0,
		"plsMethod" = 0L,
		"numThreads" = 1L,
	    "maxNComp" = 0L,
		"userEvalFunction" = object@evalFunction,
		"statistic" = 0L
	));
});

#' @rdname toCControlList-methods
setMethod("toCControlList", signature(object = "GenAlgLMEvaluator"), function(object) {
	return(list(
		"evaluatorClass" = 2L,
		"numReplications" = 0L,
	    "innerSegments" = 0L,
	    "outerSegments" = 0L,
	    "testSetSize" = 0.0,
	    "sdfact" = 0.0,
		"plsMethod" = 0L,
		"numThreads" = object@numThreads,
	    "maxNComp" = 0L,
		"userEvalFunction" = function() {NULL;},
		"statistic" = object@statisticId
	));
});

#' @rdname toCControlList-methods
setMethod("toCControlList", signature(object = "GenAlgControl"), function(object) {
	return(list(
		"populationSize" = object@populationSize,
		"numGenerations" = object@numGenerations,
		"minVariables" = object@minVariables,
		"maxVariables" = object@maxVariables,
		"elitism" = object@elitism,
		"mutationProb" = object@mutationProbability,
		"crossover" = object@crossoverId,
		"maxDuplicateEliminationTries" = object@maxDuplicateEliminationTries,
		"badSolutionThreshold" = object@badSolutionThreshold,
		"verbosity" = object@verbosity,
		"fitnessScaling" = object@fitnessScalingId
	));
});

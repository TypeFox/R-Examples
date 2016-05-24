#' Evaluator Base Class
#'
#' Virtual base class of all available evaluators
#' @aliases GenAlgEvaluator
#' @rdname GenAlgEvaluator-class
setClass("GenAlgEvaluator", representation(), contains = "VIRTUAL");

#' PLS Evaluator
#'
#' @slot numReplications The number of replications used to evaluate a variable subset.
#' @slot innerSegments The number of inner RDCV segments used in one replication.
#' @slot outerSegments The number of outer RDCV segments used in one replication.
#' @slot testSetSize The relative size of the test set (between 0 and 1).
#' @slot sdfact The factor to scale the stand. dev. of the MSEP values when selecting the optimal number
#'      of components. For the "one standard error rule", \code{sdfact} is 1.
#' @slot numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads).
#' @slot maxNComp The maximum number of components to consider in the PLS model.
#' @slot method The PLS method used to fit the PLS model (currently only SIMPLS is implemented).
#' @slot methodId The ID of the PLS method used to fit the PLS model (see C++ code for allowed values).
#' @aliases GenAlgPLSEvaluator
#' @rdname GenAlgPLSEvaluator-class
setClass("GenAlgPLSEvaluator", representation(
	numReplications = "integer",
	innerSegments = "integer",
	outerSegments = "integer",
	testSetSize = "numeric",
	sdfact = "numeric",
	numThreads = "integer",
    maxNComp = "integer",
	method = "character",
	methodId = "integer"
), validity = function(object) {
	errors <- character(0);
	MAXUINT16 <- 2^16; # unsigned 16bit integers are used (uint16_t) in the C++ code

	if(object@numThreads < 0L || object@numThreads > MAXUINT16) {
		errors <- c(errors, paste("The maximum number of threads must be greater than or equal 0 and less than", MAXUINT16));
	}

    if(object@innerSegments <= 1L || object@innerSegments > MAXUINT16) {
        errors <- c(errors, paste("The number of inner segments must be between 2 and", MAXUINT16));
    }

    if(object@maxNComp < 0L || object@maxNComp > MAXUINT16) {
        errors <- c(errors, paste("The maximum number of components must be greater than or equal 0 and less than", MAXUINT16));
    }

    if(object@testSetSize < 0 || object@testSetSize >= 1) {
        errors <- c(errors, "The test set size must be between 0 and 1.");
    }

    if(object@sdfact < 0) {
        errors <- c(errors, "The `sdfact` must be greater than or equal 0.");
    }

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
},contains = "GenAlgEvaluator");

#' User Function Evaluator
#'
#' @slot evalFunction The function that is called to evaluate the variable subset.
#' @slot sepFunction The function that calculates the standard error of prediction for the found subsets.
#' @aliases GenAlgUserEvaluator
#' @rdname GenAlgUserEvaluator-class
setClass("GenAlgUserEvaluator", representation(
	evalFunction = "function",
	sepFunction = "function"
), prototype(
	sepFunction = function(genAlg) {
		warning("Evaluator doesn't support SEP calculation -- using raw fitness");
		return(genAlg@rawFitness);
	}
), contains = "GenAlgEvaluator");

#' Fit Evaluator
#'
#' @slot numSegments The number of CV segments used in one replication.
#' @slot numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads).
#' @slot maxNComp The maximum number of components to consider in the PLS model.
#' @slot sdfact The factor to scale the stand. dev. of the MSEP values when selecting the optimal number
#'      of components. For the "one standard error rule", \code{sdfact} is 1.
#' @slot statistic The statistic used to evaluate the fitness.
#' @slot statisticId The (internal) numeric ID of the statistic.
#' @aliases GenAlgFitEvaluator
#' @rdname GenAlgFitEvaluator-class
setClass("GenAlgFitEvaluator", representation(
	numSegments = "integer",
	numThreads = "integer",
    maxNComp = "integer",
	sdfact = "numeric",
	statistic = "character",
	statisticId = "integer"
), validity = function(object) {
	errors <- character(0);
	MAXUINT16 <- 2^16; # unsigned 16bit integers are used (uint16_t) in the C++ code

	if(object@numThreads < 0L || object@numThreads > MAXUINT16) {
		errors <- c(errors, paste("The maximum number of threads must be greater than or equal 0 and less than", MAXUINT16));
	}

    if(object@numSegments <= 1L || object@numSegments > MAXUINT16) {
        errors <- c(errors, paste("The number of segments must be between 2 and", MAXUINT16));
    }

    if(object@maxNComp < 0L || object@maxNComp > MAXUINT16) {
        errors <- c(errors, paste("The maximum number of components must be greater than or equal 0 and less than", MAXUINT16));
    }

    if(object@sdfact < 0) {
        errors <- c(errors, "The `sdfact` must be greater than or equal 0.");
    }

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
},contains = "GenAlgEvaluator");


#' LM Evaluator
#'
#' @slot statistic The statistic used to evaluate the fitness.
#' @slot statisticId The (internal) numeric ID of the statistic.
#' @slot numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads).
#' @aliases GenAlgLMEvaluator
#' @rdname GenAlgLMEvaluator-class
setClass("GenAlgLMEvaluator", representation(
	statistic = "character",
	statisticId = "integer",
	numThreads = "integer"
), prototype(covariatesPC = matrix()), contains = "GenAlgEvaluator",
validity = function(object) {
	errors <- character(0);

	MAXUINT16 <- 2^16; # unsigned 16bit integers are used (uint16_t) in the C++ code

	if(object@numThreads < 0L || object@numThreads > MAXUINT16) {
		errors <- c(errors, paste("The maximum number of threads must be greater than or equal 0 and less than", MAXUINT16));
	}

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
});

#' PLS Evaluator
#'
#' Creates the object that controls the evaluation step in the genetic algorithm
#'
#' With this method the genetic algorithm uses PLS regression models to assess the prediction power of
#' variable subsets. By default, simple repeated cross-validation (srCV) is used. The optimal number
#' of PLS components is estimated using cross-validation (with \code{innerSegments} segments) on a
#' training set. The prediction power is then evaluated by fitting a PLS regression model with this optimal
#' number of components to the training set and predicting the values of a test set (of either
#' \code{testSetSize} size or \code{1 / innerSegments}, if \code{testSetSize} is not specified).
#'
#' If the parameter \code{outerSegments} is given, repeated double cross-validation is used instead.
#' There, the data set is first split into \code{outerSegments} segments and one segment is used as
#' prediction set and the other segments as test set. This is repeated for each outer segment.
#'
#' The whole procedure is repeated \code{numReplications} times to get a more reliable estimate of the
#' prediction power.
#'
#' @param numReplications The number of replications used to evaluate a variable subset (must be between 1 and 2^16)
#' @param innerSegments The number of CV segments used in one replication (must be between 2 and 2^16)
#' @param outerSegments The number of outer CV segments used in one replication (between 0 and 2^16). If this
#'      is greater than 1, repeated double cross-validation strategy (rdCV) will be used instead of
#'      simple repeated cross-validation (srCV) (see details)
#' @param testSetSize The relative size of the test set used for simple repeated CV (between 0 and 1). This parameter
#'      is ignored if outerSegments > 1 and a warning will be issued.
#' @param numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than
#'      1 or NULL means no threads)
#' @param maxNComp The maximum number of components the PLS models should consider (if not specified,
#'      the number of components is not constrained)
#' @param method The PLS method used to fit the PLS model (currently only SIMPLS is implemented)
#' @param sdfact The factor to scale the stand. dev. of the MSEP values when selecting the optimal number
#'      of components. For the "one standard error rule", \code{sdfact} is 1.
#' @return Returns an S4 object of type \code{\link{GenAlgPLSEvaluator}} to be used as argument to
#'      a call of \code{\link{genAlg}}.
#' @export
#' @family GenAlg Evaluators
#' @example examples/genAlg.R
#' @rdname GenAlgPLSEvaluator-constructor
evaluatorPLS <- function(numReplications = 30L, innerSegments = 7L, outerSegments = 1L, testSetSize = NULL,
    numThreads = NULL, maxNComp = NULL, method = c("simpls"), sdfact = 1) {
	method <- match.arg(method);

	methodId <- switch(method,
		simpls = 0L
	);

	if(is.numeric(numReplications)) {
		numReplications <- as.integer(numReplications);
	}

	if(is.numeric(innerSegments)) {
		innerSegments <- as.integer(innerSegments);
	}

	if(is.numeric(outerSegments)) {
		outerSegments <- as.integer(outerSegments);
	}

    if(outerSegments > 1 && is.numeric(testSetSize)) {
        warning("outerSegments AND testSetSize have been set. testSetSize will be ignored and rdCV with outerSegments will be used.");
    }

    if(!is.numeric(testSetSize)) {
        testSetSize <- 0;
    }

	if(missing(numThreads) || is.null(numThreads)) {
		numThreads <- 1L;
	} else if(is.numeric(numThreads)) {
		numThreads <- as.integer(numThreads);
	}

    if(missing(maxNComp) || is.null(maxNComp)) {
        maxNComp <- 0L;
    } else if(is.numeric(maxNComp)) {
        maxNComp <- as.integer(maxNComp);
    }

	return(new("GenAlgPLSEvaluator",
		numReplications = numReplications,
		innerSegments = innerSegments,
		outerSegments = outerSegments,
		testSetSize = testSetSize,
		numThreads = numThreads,
		sdfact = sdfact,
        maxNComp = maxNComp,
		method = method,
		methodId = methodId
	));
};

#' Fit Evaluator
#'
#' Creates the object that controls the evaluation step in the genetic algorithm
#'
#' The fitness of a variable subset is assessed by how well a PLS model fits the data. To estimate
#' the optimal number of components for the PLS model, cross-validation is used.
#'
#' @param numSegments The number of CV segments used to estimate the optimal number of PLS components (between 2 and 2^16).
#' @param statistic The statistic used to evaluate the fitness (BIC, AIC, adjusted R^2, or R^2).
#' @param numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than
#'      1 or NULL means no threads).
#' @param maxNComp The maximum number of components the PLS models should consider (if not specified,
#'      the number of components is not constrained)
#' @param sdfact The factor to scale the stand. dev. of the MSEP values when selecting the optimal number
#'      of components. For the "one standard error rule", \code{sdfact} is 1.
#' @return Returns an S4 object of type \code{\link{GenAlgFitEvaluator}} to be used as argument to
#'      a call of \code{\link{genAlg}}.
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorFit.R
#' @rdname GenAlgFitEvaluator-constructor
evaluatorFit <- function(numSegments = 7L, statistic = c("BIC", "AIC", "adjusted.r.squared", "r.squared"),
    numThreads = NULL, maxNComp = NULL, sdfact = 1) {
	statistic <- match.arg(statistic);

    statId = switch(statistic,
        BIC = 0L,
        AIC = 1L,
        adjusted.r.squared = 2L,
        r.squared = 3L
    );

    if(missing(numThreads) || is.null(numThreads)) {
        numThreads <- 1L;
    } else if(is.numeric(numThreads)) {
        numThreads <- as.integer(numThreads);
    }

	if(is.numeric(numSegments)) {
		numSegments <- as.integer(numSegments);
	}

	if(missing(numThreads) || is.null(numThreads)) {
		numThreads <- 1L;
	} else if(is.numeric(numThreads)) {
		numThreads <- as.integer(numThreads);
	}

    if(missing(maxNComp) || is.null(maxNComp)) {
        maxNComp <- 0L;
    } else if(is.numeric(maxNComp)) {
        maxNComp <- as.integer(maxNComp);
    }

	return(new("GenAlgFitEvaluator",
		numSegments = numSegments,
		numThreads = numThreads,
        maxNComp = maxNComp,
        sdfact = sdfact,
		statistic = statistic,
		statisticId = statId
	));
};

#' User Defined Evaluator
#'
#' Create an evaluator that uses a user defined function to evaluate the fitness
#'
#' The user specified function must take a the response vector as first and the covariates matrix as second argument.
#' The function must return a number representing the fitness of the variable subset (the higher the value the fitter the subset)
#' Additionally the user can specify a function that takes a \code{\link{GenAlg}} object and returns
#' the standard error of prediction of the found variable subsets.
#'
#' @param FUN Function used to evaluate the fitness
#' @param sepFUN Function to calculate the SEP of the variable subsets
#' @param ... Additional arguments passed to FUN and sepFUN
#' @return Returns an S4 object of type \code{\link{GenAlgUserEvaluator}}
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorUserFunction.R
#' @rdname GenAlgUserEvaluator-constructor
evaluatorUserFunction <- function(FUN, sepFUN = NULL, ...) {
	if(!is.function(FUN)) {
		stop("FUN must be of type `function`");
	};

	evalFunction <- function(y, X) {
		FUN(y, X, ...);
	};

	if(!missing(sepFUN) && is.function(sepFUN)) {
		return(new("GenAlgUserEvaluator",
			evalFunction = evalFunction,
			sepFunction = function(object, genAlg) {
				sepFUN(genAlg, ...);
			}
		));
	} else {
		return(new("GenAlgUserEvaluator",
			evalFunction = evalFunction
		));
	}
};

#' LM Evaluator
#'
#' Create an evaluator that uses a linear model to evaluate the fitness.
#'
#' Different statistics to evaluate the fitness of the variable subset can be given. If a maximum
#' absolute correlation is given the algorithm will be very slow (as the C++ implementation can not
#' be used anymore) and multithreading is not available.
#'
#' @param statistic The statistic used to evaluate the fitness
#' @param numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads)
#' @return Returns an S4 object of type \code{\link{GenAlgLMEvaluator}}
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorLM.R
#' @rdname GenAlgLMEvaluator-constructor
evaluatorLM <- function(statistic = c("BIC", "AIC", "adjusted.r.squared", "r.squared"), numThreads = NULL) {
	statistic <- match.arg(statistic);

	statId = switch(statistic,
		BIC = 0L,
		AIC = 1L,
		adjusted.r.squared = 2L,
		r.squared = 3L
	);

	if(missing(numThreads) || is.null(numThreads)) {
		numThreads <- 1L;
	} else if(is.numeric(numThreads)) {
		numThreads <- as.integer(numThreads);
	}

	return(new("GenAlgLMEvaluator",
		statistic = statistic,
		statisticId = statId,
		numThreads = numThreads
	));
};

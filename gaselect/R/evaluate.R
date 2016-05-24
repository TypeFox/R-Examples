#' Evaluate the fitness of variable subsets
#'
#' Evaluate the given variable subsets with the given Evaluator
#'
#' @param object The GenAlgEvaluator object that is used to evaluate the variables
#' @param X The data matrix used to for fitting the model
#' @param y The response vector
#' @param subsets The logical matrix where a column stands for one subset to evaluate
#' @param seed The value to seed the random number generator before evaluating
#' @param verbosity A value between 0 (no output at all) and 5 (maximum verbosity)
#' @import Rcpp
#' @useDynLib gaselect
#' @include Evaluator.R formatSegmentation.R
#' @rdname evaluate-methods
setGeneric("evaluate", function(object, X, y, subsets, seed, verbosity) { standardGeneric("evaluate"); });

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "matrix", seed = "integer", verbosity = "integer"),
function(object, X, y, subsets, seed, verbosity) {
    if(!is.logical(subsets)) {
        stop("subsets must be logical.");
    }

    if(!is.numeric(X)) {
        stop("X must be numeric.");
    }

    if(verbosity < 0L) {
        verbosity <- 0L;
    } else if(verbosity > 5L) {
        verbosity <- 5L;
    }

    if(nrow(subsets) != ncol(X)) {
        stop("The number of rows of subsets must match the number of columns of X.");
    }

    ctrlArg <- toCControlList(object);
    ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
    ctrlArg$verbosity <- verbosity;
    res <- .Call("evaluate", ctrlArg, as.matrix(X), as.matrix(y), subsets, seed, PACKAGE = "gaselect");

    res$fitness <- trueFitnessVal(object, res$fitness);
    res$segmentation <- formatSegmentation(object, res$segmentation);

    return(res);
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "logical", seed = "integer", verbosity = "integer"),
    function(object, X, y, subsets, seed, verbosity) {
    	evaluate(object, X, y, as.matrix(subsets), seed, verbosity);
    });

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "ANY", seed = "missing", verbosity = "integer"),
    function(object, X, y, subsets, seed, verbosity) {
    	evaluate(object, X, y, subsets, as.integer(sample.int(2^16, 1)), verbosity);
    });

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "ANY", seed = "integer", verbosity = "missing"),
    function(object, X, y, subsets, seed, verbosity) {
        evaluate(object, X, y, subsets, seed, 0L);
    });

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "ANY", seed = "missing", verbosity = "missing"),
    function(object, X, y, subsets, seed, verbosity) {
        evaluate(object, X, y, subsets, as.integer(sample.int(2^16, 1)), 0L);
    });


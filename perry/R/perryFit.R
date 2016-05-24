# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Resampling-based prediction error for model evaluation
#' 
#' Estimate the prediction error of a model via (repeated) \eqn{K}-fold 
#' cross-validation, (repeated) random splitting (also known as random 
#' subsampling or Monte Carlo cross-validation), or the bootstrap.  It is 
#' thereby possible to supply an object returned by a model fitting function, 
#' a model fitting function itself, or an unevaluated function call to a model 
#' fitting function.
#' 
#' (Repeated) \eqn{K}-fold cross-validation is performed in the following 
#' way.  The data are first split into \eqn{K} previously obtained blocks of 
#' approximately equal size (given by \code{folds}).  Each of the \eqn{K} data 
#' blocks is left out once to fit the model, and predictions are computed for 
#' the observations in the left-out block with \code{predictFun}.  Thus a 
#' prediction is obtained for each observation.  The response variable and the 
#' obtained predictions for all observations are then passed to the prediction 
#' loss function \code{cost} to estimate the prediction error.  For repeated 
#' \eqn{K}-fold cross-validation (as indicated by \code{splits}), this process 
#' is replicated and the estimated prediction errors from all replications are 
#' returned.
#' 
#' (Repeated) random splitting is performed similarly.  In each replication, 
#' the data are split into a training set and a test set at random.  Then the 
#' training data is used to fit the model, and predictions are computed for the 
#' test data.  Hence only the response values from the test data and the 
#' corresponding predictions are passed to the prediction loss function 
#' \code{cost}.
#' 
#' For the bootstrap estimator, each bootstrap sample is used as training data 
#' to fit the model.  The out-of-bag estimator uses the observations that do 
#' not enter the bootstrap sample as test data and computes the prediction loss 
#' function \code{cost} for those out-of-bag observations.  The 0.632 estimator 
#' is computed as a linear combination of the out-of-bag estimator and the 
#' prediction loss of the fitted values of the model computed from the full 
#' sample.
#' 
#' In any case, if the response is a vector but \code{predictFun} returns a 
#' matrix, the prediction error is computed for each column.  A typical use 
#' case for this behavior would be if \code{predictFun} returns predictions 
#' from an initial model fit and stepwise improvements thereof.
#' 
#' If \code{formula} or \code{data} are supplied, all variables required for 
#' fitting the models are added as one argument to the function call, which is 
#' the typical behavior of model fitting functions with a 
#' \code{\link[stats]{formula}} interface.  In this case, the accepted values 
#' for \code{names} depend on the method.  For the \code{function} method, a 
#' character vector of length two should supplied, with the first element 
#' specifying the argument name for the formula and the second element 
#' specifying the argument name for the data (the default is to use 
#' \code{c("formula", "data")}).  Note that names for both arguments should be 
#' supplied even if only one is actually used.  For the other methods, which do 
#' not have a \code{formula} argument, a character string specifying the 
#' argument name for the data should be supplied (the default is to use 
#' \code{"data"}).  
#' 
#' If \code{x} is supplied, on the other hand, the predictor matrix and the 
#' response are added as separate arguments to the function call.  In this 
#' case, \code{names} should be a character vector of length two, with the 
#' first element specifying the argument name for the predictor matrix and the 
#' second element specifying the argument name for the response (the default is 
#' to use \code{c("x", "y")}).  It should be noted that the \code{formula} or 
#' \code{data} arguments take precedence over \code{x}.
#' 
#' @aliases print.perry
#' 
#' @param object  the fitted model for which to estimate the prediction error, 
#' a function for fitting a model, or an unevaluated function call for fitting 
#' a model (see \code{\link{call}} for the latter).  In the case of a fitted 
#' model, the object is required to contain a component \code{call} that stores 
#' the function call used to fit the model, which is typically the case for 
#' objects returned by model fitting functions.
#' @param formula  a \code{\link[stats]{formula}} describing the model.
#' @param data  a data frame containing the variables required for fitting the 
#' models.  This is typically used if the model in the function call is 
#' described by a \code{\link[stats]{formula}}.
#' @param x  a numeric matrix containing the predictor variables.  This is 
#' typically used if the function call for fitting the models requires the 
#' predictor matrix and the response to be supplied as separate arguments.
#' @param y  a numeric vector or matrix containing the response.
#' @param args  a list of additional arguments to be passed to the model 
#' fitting function.
#' @param splits  an object of class \code{"cvFolds"} (as returned by 
#' \code{\link{cvFolds}}) or a control object of class \code{"foldControl"} 
#' (see \code{\link{foldControl}}) defining the folds of the data for 
#' (repeated) \eqn{K}-fold cross-validation, an object of class 
#' \code{"randomSplits"} (as returned by \code{\link{randomSplits}}) or a 
#' control object of class \code{"splitControl"} (see 
#' \code{\link{splitControl}}) defining random data splits, or an object of 
#' class \code{"bootSamples"} (as returned by \code{\link{bootSamples}}) or a 
#' control object of class \code{"bootControl"} (see \code{\link{bootControl}}) 
#' defining bootstrap samples.
#' @param predictFun  a function to compute predictions for the test data.  It 
#' should expect the fitted model to be passed as the first argument and the test 
#' data as the second argument, and must return either a vector or a matrix 
#' containing the predicted values.  The default is to use the 
#' \code{\link[stats]{predict}} method of the fitted model.
#' @param predictArgs  a list of additional arguments to be passed to  
#' \code{predictFun}.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' the observed values of the response to be passed as the first argument and 
#' the predicted values as the second argument, and must return either a 
#' non-negative scalar value, or a list with the first component containing 
#' the prediction error and the second component containing the standard 
#' error.  The default is to use the root mean squared prediction error 
#' (see \code{\link{cost}}).
#' @param costArgs  a list of additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' @param names  an optional character vector giving names for the arguments 
#' containing the data to be used in the function call (see \dQuote{Details}).
#' @param envir  the \code{\link{environment}} in which to evaluate the 
#' function call for fitting the models (see \code{\link{eval}}).
#' @param ncores  a positive integer giving the number of processor cores to be 
#' used for parallel computing (the default is 1 for no parallelization).  If 
#' this is set to \code{NA}, all available processor cores are used.
#' @param cl  a \pkg{parallel} cluster for parallel computing as generated by 
#' \code{\link[parallel]{makeCluster}}.  If supplied, this is preferred over 
#' \code{ncores}.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).  Note that also in case of parallel computing, 
#' resampling is performed on the manager process rather than the worker 
#' processes. On the parallel worker processes, random number streams are 
#' used and the seed is set via \code{\link{clusterSetRNGStream}} for 
#' reproducibility in case the model fitting function involves randomness.
#' @param \dots  additional arguments to be passed down.
#' 
#' @returnClass perry
#' @returnItem pe  a numeric vector containing the respective estimated 
#' prediction errors.  In case of more than one replication, those are average 
#' values over all replications.
#' @returnItem se  a numeric vector containing the respective estimated 
#' standard errors of the prediction loss.
#' @returnItem reps  a numeric matrix in which each column contains the 
#' respective estimated prediction errors from all replications.  This is 
#' only returned in case of more than one replication.
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error.
#' @returnItem y  the response.
#' @returnItem yHat  a list containing the predicted values from all 
#' replications.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perrySelect}}, \code{\link{perryTuning}}, 
#' \code{\link{cvFolds}}, \code{\link{randomSplits}}, 
#' \code{\link{bootSamples}}, \code{\link{cost}}
#' 
#' @example inst/doc/examples/example-perryFit.R
#' 
#' @keywords utilities
#' 
#' @import parallel
#' @export

perryFit <- function(object, ...) UseMethod("perryFit")


#' @rdname perryFit
#' @method perryFit default
#' @export

perryFit.default <- function(object, data = NULL, x = NULL, y, 
        splits = foldControl(), predictFun = predict, predictArgs = list(), 
        cost = rmspe, costArgs = list(), names = NULL, envir = parent.frame(), 
        ncores = 1, cl = NULL, seed = NULL, ...) {
    ## extract function call for model fit
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    call <- object$call
    if(is.null(call)) stop("function call for model fitting not available")
    ## for the 0.632 bootstrap estimator, compute the apparent error
    if(inherits(splits, c("bootControl", "bootSamples")) && splits$type == "0.632") {
        # predict the response for all observations
        if(is.null(data)) 
            splits$yHat <- doCall(predictFun, object, x, args=predictArgs)
        else splits$yHat <- doCall(predictFun, object, data, args=predictArgs)
    }
    ## call method for unevaluated function calls
    out <- perryFit(call, data=data, x=x, y=y, splits=splits, 
        predictFun=predictFun, predictArgs=predictArgs, cost=cost, 
        costArgs=costArgs, names=names, envir=envir, ncores=ncores, 
        cl=cl, seed=seed, ...)
    out$call <- matchedCall
    out
}


#' @rdname perryFit
#' @method perryFit function
#' @export

perryFit.function <- function(object, formula, data = NULL, x = NULL, y, 
        args = list(), splits = foldControl(), predictFun = predict, 
        predictArgs = list(), cost = rmspe, costArgs = list(), names = NULL, 
        envir = parent.frame(), ncores = 1, cl = NULL, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    call <- as.call(c(object, args))  # set up unevaluated function call
    haveFormula <- !missing(formula)
    if(haveFormula || !missing(data)) {
        if(is.null(names)) names <- c("formula", "data")
        if(haveFormula) call[[names[1]]] <- formula
        names <- names[-1]
        mf <- match.call(expand.dots = FALSE)
        m <- match(c("formula", "data"), names(mf), 0)
        mf <- mf[c(1, m)]
        mf$drop.unused.levels <- TRUE
        mf[[1]] <- as.name("model.frame")
        data <- eval(mf, envir)
        if(is.empty.model(attr(data, "terms"))) stop("empty model")
        y <- model.response(data)  # extract response from model frame
    }
    ## call method for unevaluated function calls
    out <- perryFit(call, data=data, x=x, y=y, splits=splits, 
        predictFun=predictFun, predictArgs=predictArgs, cost=cost, 
        costArgs=costArgs, names=names, envir=envir, ncores=ncores, 
        cl=cl, seed=seed, ...)
    out$call <- matchedCall
    out
}


#' @rdname perryFit
#' @method perryFit call
#' @export

perryFit.call <- function(object, data = NULL, x = NULL, y, 
        splits = foldControl(), predictFun = predict, predictArgs = list(), 
        cost = rmspe, costArgs = list(), names = NULL, envir = parent.frame(), 
        ncores = 1, cl = NULL, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("perryFit")
    n <- nobs(y)
    if(is.null(data)) {
        sx <- "x"
        nx <- nobs(x)
    } else {
        sx <- "data"
        nx <- nobs(data)
    }
    if(!isTRUE(n == nx)) stop(sprintf("'%s' must have %d observations", sx, nx))
    if(!is.null(seed)) set.seed(seed)
    ## compute data splits
    if(hasMethod("perrySplits", class(splits))) splits <- perrySplits(n, splits)
    ## compute fitted values from the model using all observations for the 
    ## 0.632 bootstrap estimator (if they are not yet available)
    if(inherits(splits, "bootSamples") && splits$type == "0.632") {
        # check if the apparent error has already been computed
        if(!hasComponent(splits, "yHat")) {
            # fit the model from all observations
            if(is.null(data)) {
                if(is.null(names)) names <- c("x", "y")
                object[[names[1]]] <- x
                object[[names[2]]] <- y
            } else {
                if(is.null(names)) names <- "data"
                object[[names]] <- data
            }
            fit <- eval(object, envir)
            # predict the response for all observations
            if(is.null(data)) 
                splits$yHat <- doCall(predictFun, fit, x, args=predictArgs)
            else splits$yHat <- doCall(predictFun, fit, data, args=predictArgs)
        }
    }
    # set up parallel computing if requested
    R <- splits$R
    haveCl <- inherits(cl, "cluster")
    if(haveCl) haveNcores <- FALSE
    else {
        if(is.na(ncores)) ncores <- detectCores()  # use all available cores
        if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
            ncores <- 1  # use default value
            warning("invalid value of 'ncores'; using default value")
        } else ncores <- as.integer(ncores)
        if(inherits(splits, "cvFolds") && R == 1) {
            ncores <- min(ncores, splits$K)
        } else ncores <- min(ncores, R)
        haveNcores <- ncores > 1
    }
    # check whether parallel computing should be used
    useParallel <- haveNcores || haveCl
    # set up multicore or snow cluster if not supplied
    if(haveNcores) {
        if(.Platform$OS.type == "windows") {
            cl <- makePSOCKcluster(rep.int("localhost", ncores))
        } else cl <- makeForkCluster(ncores)
        on.exit(stopCluster(cl))
    }
    if(useParallel) {
        # set seed of the random number stream
        if(!is.null(seed)) clusterSetRNGStream(cl, iseed=seed)
        else if(haveNcores) clusterSetRNGStream(cl)
    }
    ## call workhorse function to compute the predictions
    yHat <- perryPredictions(object, data, x, y, splits=splits, 
        predictFun=predictFun, predictArgs=predictArgs, names=names, 
        envir=envir, cl=cl)
    ## call workhorse function to estimate the prediction loss
    pe <- perryCost(splits, y, yHat, cost=cost, costArgs=costArgs)
    ## construct return object
    pe <- c(pe, list(splits=splits, y=y, yHat=yHat, call=matchedCall))
    class(pe) <- "perry"
    pe
}

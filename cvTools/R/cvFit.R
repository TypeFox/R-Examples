# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Cross-validation for model evaluation
#' 
#' Estimate the prediction error of a model via (repeated) \eqn{K}-fold 
#' cross-validation.  It is thereby possible to supply an object returned by a 
#' model fitting function, a model fitting function itself, or an unevaluated 
#' function call to a model fitting function.
#' 
#' (Repeated) \eqn{K}-fold cross-validation is performed in the following 
#' way.  The data are first split into \eqn{K} previously obtained blocks of 
#' approximately equal size.  Each of the \eqn{K} data blocks is left out once 
#' to fit the model, and predictions are computed for the observations in the 
#' left-out block with the \code{\link[stats]{predict}} method of the fitted 
#' model.  Thus a prediction is obtained for each observation.
#' 
#' The response variable and the obtained predictions for all observations are 
#' then passed to the prediction loss function \code{cost} to estimate the 
#' prediction error.  For repeated cross-validation, this process is replicated 
#' and the estimated prediction errors from all replications as well as their 
#' average are included in the returned object.
#' 
#' Furthermore, if the response is a vector but the 
#' \code{\link[stats]{predict}} method of the fitted models returns a matrix, 
#' the prediction error is computed for each column.  A typical use case for 
#' this behavior would be if the \code{\link[stats]{predict}} method returns 
#' predictions from an initial model fit and stepwise improvements thereof.
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
#' @aliases print.cv
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
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' the observed values of the response to be passed as the first argument and 
#' the predicted values as the second argument, and must return either a 
#' non-negative scalar value, or a list with the first component containing 
#' the prediction error and the second component containing the standard 
#' error.  The default is to use the root mean squared prediction error 
#' (see \code{\link{cost}}).
#' @param K  an integer giving the number of groups into which the data should 
#' be split (the default is five).  Keep in mind that this should be chosen 
#' such that all groups are of approximately equal size.  Setting \code{K} 
#' equal to \code{n} yields leave-one-out cross-validation.
#' @param R  an integer giving the number of replications for repeated 
#' \eqn{K}-fold cross-validation.  This is ignored for for leave-one-out 
#' cross-validation and other non-random splits of the data.
#' @param foldType  a character string specifying the type of folds to be 
#' generated.  Possible values are \code{"random"} (the default), 
#' \code{"consecutive"} or \code{"interleaved"}.
#' @param folds  an object of class \code{"cvFolds"} giving the folds of the 
#' data for cross-validation (as returned by \code{\link{cvFolds}}).  If 
#' supplied, this is preferred over \code{K} and \code{R}.
#' @param names  an optional character vector giving names for the arguments 
#' containing the data to be used in the function call (see \dQuote{Details}).
#' @param predictArgs  a list of additional arguments to be passed to the 
#' \code{\link[stats]{predict}} method of the fitted models.
#' @param costArgs  a list of additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' @param envir  the \code{\link{environment}} in which to evaluate the 
#' function call for fitting the models (see \code{\link{eval}}).
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed down.
#' 
#' @returnClass cv
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem cv  a numeric vector containing the respective estimated 
#' prediction errors.  For repeated cross-validation, those are average values 
#' over all replications.
#' @returnItem se  a numeric vector containing the respective estimated 
#' standard errors of the prediction loss.
#' @returnItem reps  a numeric matrix in which each column contains the 
#' respective estimated prediction errors from all replications.  This is 
#' only returned for repeated cross-validation.
#' @returnItem seed  the seed of the random number generator before 
#' cross-validation was performed.
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvTool}}, \code{\link{cvSelect}}, 
#' \code{\link{cvTuning}}, \code{\link{cvFolds}}, \code{\link{cost}}
#' 
#' @example inst/doc/examples/example-cvFit.R
#' 
#' @keywords utilities
#' 
#' @export

cvFit <- function(object, ...) UseMethod("cvFit")


#' @rdname cvFit
#' @method cvFit default
#' @export

cvFit.default <- function(object, data = NULL, x = NULL, y, cost = rmspe, 
        K = 5, R = 1, foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, names = NULL, predictArgs = list(), costArgs = list(), 
        envir = parent.frame(), seed = NULL, ...) {
    ## extract function call for model fit
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("cvFit")
    call <- object$call
    if(is.null(call)) stop("function call for model fitting not available")
    ## call method for unevaluated function calls
    out <- cvFit(call, data=data, x=x, y=y, cost=cost, K=K, R=R, 
        foldType=foldType, folds=folds, names=names, predictArgs=predictArgs, 
        costArgs=costArgs, envir=envir, seed=seed)
    out$call <- matchedCall
    out
}


#' @rdname cvFit
#' @method cvFit function
#' @export

cvFit.function <- function(object, formula, data = NULL, x = NULL, y, 
        args = list(), cost = rmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), folds = NULL, 
        names = NULL, predictArgs = list(), costArgs = list(), 
        envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("cvFit")
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
    out <- cvFit(call, data=data, x=x, y=y, cost=cost, K=K, R=R, 
        foldType=foldType, folds=folds, names=names, predictArgs=predictArgs, 
        costArgs=costArgs, envir=envir, seed=seed)
    out$call <- matchedCall
    out
}


#' @rdname cvFit
#' @method cvFit call
#' @export

cvFit.call <- function(object, data = NULL, x = NULL, y, cost = rmspe, 
        K = 5, R = 1, foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, names = NULL, predictArgs = list(), costArgs = list(), 
        envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("cvFit")
    n <- nobs(y)
    if(is.null(data)) {
        sx <- "x"
        nx <- nobs(x)
    } else {
        sx <- "data"
        nx <- nobs(data)
    }
    if(!isTRUE(n == nx)) stop(sprintf("'%s' must have %d observations", sx, nx))
    # make sure that .Random.seed exists if no seed is supplied
    if(is.null(seed)) {
        if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
        seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    } else set.seed(seed)
    ## compute data blocks as in 'cv.lars'
    if(is.null(folds)) folds <- cvFolds(n, K, R, type=foldType)
    R <- folds$R
    ## call workhorse function to perform cross-validation
    cv <- cvTool(object, data, x, y, cost=cost, folds=folds, names=names, 
        predictArgs=predictArgs, costArgs=costArgs, envir=envir)
    ## compute average results in case of repeated CV
    if(R > 1) {
        reps <- cv
        cv <- apply(reps, 2, mean)
        se <- apply(reps, 2, sd)
    } else {
        if(is.list(cv)) {
            se <- cv[[2]]
            cv <- cv[[1]]
        } else {
            cv <- drop(cv)
            if(is.null(names(cv))) {
                # drop() removes column name of 1x1 matrix
                names(cv) <- defaultCvNames(length(cv))
            }
            se <- rep.int(NA, length(cv))
            names(se) <- names(cv)
        }
    }
    ## construct return object
    out <- list(n=folds$n, K=folds$K, R=R, cv=cv, se=se)
    if(R > 1) out$reps <- reps
    out$seed <- seed
    out$call <- matchedCall
    class(out) <- "cv"
    out
}

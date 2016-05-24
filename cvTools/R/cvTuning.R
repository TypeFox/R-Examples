# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Cross-validation for tuning parameter selection
#' 
#' Select tuning parameters of a model by estimating the respective prediction 
#' errors via (repeated) \eqn{K}-fold cross-validation.  It is thereby possible 
#' to supply a model fitting function or an unevaluated function call to a 
#' model fitting function.
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
#' supplied even if only one is actually used.  For the \code{call} method, 
#' which does not have a \code{formula} argument, a character string specifying 
#' the argument name for the data should be supplied (the default is to use 
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
#' @aliases print.cvTuning
#' 
#' @param object  a function or an unevaluated function call for fitting 
#' a model (see \code{\link{call}} for the latter).
#' @param formula  a \code{\link[stats]{formula}} describing the model.
#' @param data  a data frame containing the variables required for fitting the 
#' models.  This is typically used if the model in the function call is 
#' described by a \code{\link[stats]{formula}}.
#' @param x  a numeric matrix containing the predictor variables.  This is 
#' typically used if the function call for fitting the models requires the 
#' predictor matrix and the response to be supplied as separate arguments.
#' @param y  a numeric vector or matrix containing the response.
#' @param tuning  a list of arguments giving the tuning parameter values to be 
#' evaluated.  The names of the list components should thereby correspond to 
#' the argument names of the tuning parameters.  For each tuning parameter, a 
#' vector of values can be supplied.  Cross-validation is then applied over the 
#' grid of all possible combinations of tuning parameter values.
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
#' @param selectBest  a character string specifying a criterion for selecting 
#' the best model.  Possible values are \code{"min"} (the default) or 
#' \code{"hastie"}.  The former selects the model with the smallest prediction 
#' error.  The latter is useful for models with a tuning parameter controlling 
#' the complexity of the model (e.g., penalized regression).  It selects the 
#' most parsimonious model whose prediction error is no larger than 
#' \code{seFactor} standard errors above the prediction error of the best 
#' overall model.  Note that the models are thereby assumed to be ordered 
#' from the most parsimonious one to the most complex one.  In particular 
#' a one-standard-error rule is frequently applied.
#' @param seFactor  a numeric value giving a multiplication factor of the 
#' standard error for the selection of the best model.  This is ignored if 
#' \code{selectBest} is \code{"min"}.
#' @param envir  the \code{\link{environment}} in which to evaluate the 
#' function call for fitting the models (see \code{\link{eval}}).
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed down.
#' 
#' @return 
#' If \code{tuning} is an empty list, \code{\link{cvFit}} is called to return 
#' an object of class \code{"cv"}.
#' 
#' Otherwise an object of class \code{"cvTuning"} (which inherits from class 
#' \code{"cvSelect"}) with the following components is returned:
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem tuning  a data frame containing the grid of tuning parameter 
#' values for which the prediction error was estimated.
#' @returnItem best  an integer vector giving the indices of the optimal 
#' combinations of tuning parameters.
#' @returnItem cv  a data frame containing the estimated prediction errors for 
#' all combinations of tuning parameter values.  For repeated cross-validation, 
#' those are average values over all replications.
#' @returnItem se  a data frame containing the estimated standard errors of the 
#' prediction loss for all combinations of tuning parameter values.
#' @returnItem selectBest  a character string specifying the criterion used for 
#' selecting the best model.
#' @returnItem seFactor  a numeric value giving the multiplication factor of 
#' the standard error used for the selection of the best model.
#' @returnItem reps  a data frame containing the estimated prediction errors 
#' from all replications for all combinations of tuning parameter values.  This 
#' is only returned for repeated cross-validation.
#' @returnItem seed  the seed of the random number generator before 
#' cross-validation was performed.
#' @returnItem call  the matched function call.
#' 
#' @note The same cross-validation folds are used for all combinations of 
#' tuning parameter values for maximum comparability.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Hastie, T., Tibshirani, R. and Friedman, J. (2009) \emph{The Elements of 
#' Statistical Learning: Data Mining, Inference, and Prediction}.  Springer, 
#' 2nd edition.
#' 
#' @seealso \code{\link{cvTool}}, \code{\link{cvFit}}, \code{\link{cvSelect}}, 
#' \code{\link{cvFolds}}, \code{\link{cost}}
#' 
#' @example inst/doc/examples/example-cvTuning.R
#' 
#' @keywords utilities
#' 
#' @export

cvTuning <- function(object, ...) UseMethod("cvTuning")


#' @rdname cvTuning
#' @method cvTuning function
#' @export

cvTuning.function <- function(object, formula, data = NULL, x = NULL, y, 
        tuning = list(), args = list(), cost = rmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), folds = NULL, 
        names = NULL, predictArgs = list(), costArgs = list(), 
        selectBest = c("min", "hastie"), seFactor = 1, 
        envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("cvTuning")
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
    out <- cvTuning(call, data=data, x=x, y=y, tuning=tuning, cost=cost, 
        K=K, R=R, foldType=foldType, folds=folds, names=names, 
        predictArgs=predictArgs, costArgs=costArgs, selectBest=selectBest, 
        seFactor=seFactor, envir=envir, seed=seed)
    out$call <- matchedCall
    out
}


#' @rdname cvTuning
#' @method cvTuning call
#' @export

cvTuning.call <- function(object, data = NULL, x = NULL, y, 
        tuning = list(), cost = rmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, names = NULL, predictArgs = list(), 
        costArgs = list(), selectBest = c("min", "hastie"), 
        seFactor = 1, envir = parent.frame(), seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    matchedCall[[1]] <- as.name("cvTuning")
    n <- nobs(y)
    if(is.null(data)) {
        sx <- "x"
        nx <- nobs(x)
    } else {
        sx <- "data"
        nx <- nobs(data)
    }
    if(!isTRUE(n == nx)) stop(sprintf("'%s' must have %d observations", sx, nx))
    selectBest <- match.arg(selectBest)
    # create all combinations of tuning parameters
    tuning <- do.call("expand.grid", tuning)
    nTuning <- nrow(tuning)
    pTuning <- ncol(tuning)
    if(nTuning == 0 || pTuning == 0) {
        # use function cvFit() if no tuning parameters are supplied
        out <- cvFit(object, data, x, y, cost=cost, K=K, R=R, 
            foldType=foldType, folds=folds, names=names, 
            predictArgs=predictArgs, costArgs=costArgs, 
            envir=envir, seed=seed)
        return(out)
    }
    # make sure that .Random.seed exists if no seed is supplied
    if(is.null(seed)) {
        if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) runif(1)
        seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    } else set.seed(seed)
    ## compute data blocks as in 'cv.lars'
    if(is.null(folds)) folds <- cvFolds(n, K, R, type=foldType)
    R <- folds$R
    ## perform cross-validation for each combination of tuning parameters
    tuningNames <- names(tuning)
    if(pTuning == 1) {
        cv <- lapply(tuning[, tuningNames], 
            function(t) {
                object[[tuningNames]] <- t
                cvTool(object, data, x, y, cost=cost, folds=folds, names=names, 
                    predictArgs=predictArgs, costArgs=costArgs, envir=envir)
            })
    } else {
        cv <- lapply(seq_len(nTuning), 
            function(i) {
                for(j in seq_len(pTuning)) {
                    object[[tuningNames[j]]] <- tuning[i, j]
                }
                cvTool(object, data, x, y, cost=cost, folds=folds, names=names, 
                    predictArgs=predictArgs, costArgs=costArgs, envir=envir)
            })
    }
    if(R == 1) {
        haveList <- is.list(cv[[1]])
        if(haveList) {
            se <- lapply(cv, function(x) x[[2]])
            cv <- lapply(cv, function(x) x[[1]])
        }
    }
    cv <- do.call("rbind", cv)
    if(R == 1) {
        if(haveList) {
            se <- do.call("rbind", se)
        } else se <- matrix(NA, nrow(cv), ncol(cv), dimnames=dimnames(cv))
        se <- data.frame(Fit=seq_len(nTuning), se, row.names=NULL)
    }
    cv <- data.frame(Fit=rep(seq_len(nTuning), each=R), cv, row.names=NULL)
    ## compute average results in case of repeated CV
    if(R > 1) {
        reps <- cv
        cv <- aggregate(reps[, -1, drop=FALSE], reps[, 1, drop=FALSE], mean)
        se <- aggregate(reps[, -1, drop=FALSE], reps[, 1, drop=FALSE], sd)
    }
    ## find optimal values of tuning parameters
    if(selectBest == "min") {
        seFactor <- NA
        best <- sapply(cv[, -1, drop=FALSE], selectMin)
    } else {
        seFactor <- rep(seFactor, length.out=1)
        best <- sapply(names(cv)[-1], 
            function(j) selectHastie(cv[, j], se[, j], seFactor=seFactor))
    }
    ## construct return object
    out <- list(n=folds$n, K=folds$K, R=R, tuning=tuning, best=best, 
        cv=cv, se=se, selectBest=selectBest, seFactor=seFactor)
    if(R > 1) out$reps <- reps
    out$seed <- seed
    out$call <- matchedCall
    class(out) <- c("cvTuning", "cvSelect")
    out
}

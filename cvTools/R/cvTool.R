# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Low-level function for cross-validation
#' 
#' Basic function to estimate the prediction error of a model via (repeated) 
#' \eqn{K}-fold cross-validation.  The model is thereby specified by an 
#' unevaluated function call to a model fitting function.
#' 
#' (Repeated) \eqn{K}-fold cross-validation is performed in the following 
#' way.  The data are first split into \eqn{K} previously obtained blocks of 
#' approximately equal size (given by \code{folds}).  Each of the \eqn{K} data 
#' blocks is left out once to fit the model, and predictions are computed for 
#' the observations in the left-out block with the \code{\link[stats]{predict}} 
#' method of the fitted model.  Thus a prediction is obtained for each 
#' observation.
#' 
#' The response variable and the obtained predictions for all observations are 
#' then passed to the prediction loss function \code{cost} to estimate the 
#' prediction error.  For repeated cross-validation (as indicated by 
#' \code{folds}), this process is replicated and the estimated prediction 
#' errors from all replications are returned.
#' 
#' Furthermore, if the response is a vector but the 
#' \code{\link[stats]{predict}} method of the fitted models returns a matrix, 
#' the prediction error is computed for each column.  A typical use case for 
#' this behavior would be if the \code{\link[stats]{predict}} method returns 
#' predictions from an initial model fit and stepwise improvements thereof.
#' 
#' If \code{data} is supplied, all variables required for fitting the models 
#' are added as one argument to the function call, which is the typical 
#' behavior of model fitting functions with a \code{\link[stats]{formula}} 
#' interface.  In this case, a character string specifying the argument name 
#' can be passed via \code{names} (the default is to use \code{"data"}).  
#' 
#' If \code{x} is supplied, on the other hand, the predictor matrix and the 
#' response are added as separate arguments to the function call.  In this 
#' case, \code{names} should be a character vector of length two, with the 
#' first element specifying the argument name for the predictor matrix and the 
#' second element specifying the argument name for the response (the default is 
#' to use \code{c("x", "y")}).  It should be noted that \code{data} takes 
#' precedence over \code{x} if both are supplied.
#' 
#' @param call  an unevaluated function call for fitting a model (see 
#' \code{\link{call}}).
#' @param data  a data frame containing the variables required for fitting the 
#' models.  This is typically used if the model in the function call is 
#' described by a \code{\link[stats]{formula}}.
#' @param x  a numeric matrix containing the predictor variables.  This is 
#' typically used if the function call for fitting the models requires the 
#' predictor matrix and the response to be supplied as separate arguments.
#' @param y  a numeric vector or matrix containing the response.
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' the observed values of the response to be passed as the first argument and 
#' the predicted values as the second argument, and must return either a 
#' non-negative scalar value, or a list with the first component containing 
#' the prediction error and the second component containing the standard 
#' error.  The default is to use the root mean squared prediction error 
#' (see \code{\link{cost}}).
#' @param folds  an object of class \code{"cvFolds"} giving the folds of the 
#' data for cross-validation (as returned by \code{\link{cvFolds}}).
#' @param names  an optional character vector giving names for the arguments 
#' containing the data to be used in the function call (see \dQuote{Details}).
#' @param predictArgs  a list of additional arguments to be passed to the 
#' \code{\link[stats]{predict}} method of the fitted models.
#' @param costArgs  a list of additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' @param envir  the \code{\link{environment}} in which to evaluate the 
#' function call for fitting the models (see \code{\link{eval}}).
#' 
#' @return If only one replication is requested and the prediction loss 
#' function \code{cost} also returns the standard error, a list is returned, 
#' with the first component containing the estimated prediction errors and the 
#' second component the corresponding estimated standard errors.
#' 
#' Otherwise the return value is a numeric matrix in which each column contains 
#' the respective estimated prediction errors from all replications.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvTuning}}, \code{\link{cvFolds}}, 
#' \code{\link{cost}}
#' 
#' @example inst/doc/examples/example-cvTool.R
#' 
#' @keywords utilities
#' 
#' @export

cvTool <- function(call, data = NULL, x = NULL, y, cost = rmspe, folds, 
        names = NULL, predictArgs = list(), costArgs = list(),  
        envir = parent.frame()) {
    # initializations
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        cvExpr <- expression(
            lapply(subsetList, cvXY, call=call, x=x, y=y, 
                names=names, predictArgs=predictArgs, envir=envir)
        )
    } else {
        if(is.null(names)) names <- "data"
        cvExpr <- expression(
            lapply(subsetList, cvData, call=call, data=data, 
                names=names, predictArgs=predictArgs, envir=envir)
        )
    }
    R <- folds$R
    # perform (repeated) cross-validation
    if(R == 1) {
        subsetList <- getSubsetList(folds)
        yHat <- eval(cvExpr)  # perform cross-validation
        # instead of collecting the results from the folds in the original 
        # order of the observations, the response is re-ordered accordingly
        yHat <- combineData(yHat)  # combine predictions from the folds
        y <- dataSubset(y, unlist(subsetList))  # re-order response
        # compute cost function for predicted values
        if(is.null(dim(y)) && !is.null(dim(yHat))) {
            cv <- apply(yHat, 2, 
                function(yHat) doCall(cost, y, yHat, args=costArgs))
            if(is.list(cv)) {
                cv <- list(sapply(cv, function(x) x[[1]]), 
                    sapply(cv, function(x) x[[2]]))
            }
        } else cv <- doCall(cost, y, yHat, args=costArgs)
    } else {
        cv <- sapply(seq_len(R), 
            function(r) {
                subsetList <- getSubsetList(folds, r)
                yHat <- eval(cvExpr)  # perform cross-validation
                # instead of collecting the results from the folds in the original 
                # order of the observations, the response is re-ordered accordingly
                yHat <- combineData(yHat)  # combine predictions from the folds
                y <- dataSubset(y, unlist(subsetList))  # re-order response
                # compute cost function for predicted values
                if(is.null(dim(y)) && !is.null(dim(yHat))) {
                    tmp <- apply(yHat, 2, 
                        function(yHat) doCall(cost, y, yHat, args=costArgs))
                    if(is.list(tmp)) tmp <- sapply(tmp, function(x) x[[1]])
                } else {
                    tmp <- doCall(cost, y, yHat, args=costArgs)
                    if(is.list(tmp)) tmp <- tmp[[1]]
                }
                tmp
            })
    }
    if(is.list(cv)) {
        cv <- lapply(cv, 
            function(x) {
                if(is.null(names(x))) names(x) <- defaultCvNames(length(x))
                x
            })
    } else {
        cv <- if(is.null(dim(cv)) && R > 1) as.matrix(cv) else t(cv)
        if(is.null(colnames(cv))) colnames(cv) <- defaultCvNames(ncol(cv))
    }
    cv
}

cvXY <- function(i, call, x, y, names, predictArgs, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, -i)
    call[[names[2]]] <- dataSubset(y, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(x, i), args=predictArgs)
}

cvData <- function(i, call, data, names, predictArgs, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predict, fit, dataSubset(data, i), args=predictArgs)
}

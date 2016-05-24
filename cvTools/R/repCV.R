# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Cross-validation for linear models
#' 
#' Estimate the prediction error of a linear model via (repeated) \eqn{K}-fold 
#' cross-validation.  Cross-validation functions are available for least 
#' squares fits computed with \code{\link[stats]{lm}} as well as for the 
#' following robust alternatives: MM-type models computed with 
#' \code{\link[robustbase]{lmrob}} and least trimmed squares fits computed with 
#' \code{\link[robustbase]{ltsReg}}.
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
#' @aliases cvExamples
#' 
#' @param object  an object returned from a model fitting function.  Methods 
#' are implemented for objects of class \code{"lm"} computed with 
#' \code{\link[stats]{lm}}, objects of class \code{"lmrob"} computed with 
#' \code{\link[robustbase]{lmrob}}, and object of class \code{"lts"} computed 
#' with \code{\link[robustbase]{ltsReg}}. 
#' @param cost  a cost function measuring prediction loss.  It should expect 
#' the observed values of the response to be passed as the first argument and 
#' the predicted values as the second argument, and must return either a 
#' non-negative scalar value, or a list with the first component containing 
#' the prediction error and the second component containing the standard 
#' error.  The default is to use the root mean squared prediction error 
#' for the \code{"lm"} method and the root trimmed mean squared prediction 
#' error for the \code{"lmrob"} and \code{"lts"} methods (see 
#' \code{\link{cost}}).
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
#' @param fit  a character string specifying for which fit to estimate the 
#' prediction error.  Possible values are \code{"reweighted"} (the default) for 
#' the prediction error of the reweighted fit, \code{"raw"} for the prediction 
#' error of the raw fit, or \code{"both"} for the prediction error of both 
#' fits.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).
#' @param \dots  additional arguments to be passed to the prediction loss 
#' function \code{cost}.
#' 
#' @returnClass cv
#' @returnItem n  an integer giving the number of observations.
#' @returnItem K  an integer giving the number of folds.
#' @returnItem R  an integer giving the number of replications.
#' @returnItem cv  a numeric vector containing the estimated prediction 
#' errors.  For the \code{"lm"} and \code{"lmrob"} methods, this is a single 
#' numeric value.  For the \code{"lts"} method, this contains one value for 
#' each of the requested fits.  In the case of repeated cross-validation, those 
#' are average values over all replications.
#' @returnItem se  a numeric vector containing the estimated standard 
#' errors of the prediction loss.  For the \code{"lm"} and \code{"lmrob"} 
#' methods, this is a single numeric value.  For the \code{"lts"} method, this 
#' contains one value for each of the requested fits.
#' @returnItem reps  a numeric matrix containing the estimated prediction 
#' errors from all replications.  For the \code{"lm"} and \code{"lmrob"} 
#' methods, this is a matrix with one column.  For the \code{"lts"} method, 
#' this contains one column for each of the requested fits.  However, this is 
#' only returned for repeated cross-validation.
#' @returnItem seed  the seed of the random number generator before 
#' cross-validation was performed.
#' @returnItem call  the matched function call.
#' 
#' @note The \code{repCV} methods are simple wrapper functions that extract the 
#' data from the fitted model and call \code{\link{cvFit}} to perform 
#' cross-validation.  In addition, \code{cvLm}, \code{cvLmrob} and \code{cvLts} 
#' are aliases for the respective methods.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{cvFit}}, \code{\link{cvFolds}}, \code{\link{cost}}, 
#' \code{\link[stats]{lm}}, \code{\link[robustbase]{lmrob}}, 
#' \code{\link[robustbase]{ltsReg}}
#' 
#' @example inst/doc/examples/example-repCV.R
#' 
#' @keywords utilities
#' 
#' @export repCV

repCV <- function(object, ...) UseMethod("repCV")


## LS regression 
#' @rdname repCV
#' @method repCV lm
#' @export

repCV.lm <- function(object, cost=rmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    # retrieve data from model fit
    if(is.null(data <- object$model)) {
        haveDataArgument <- !is.null(object$call$data)
        if(haveDataArgument) {
            # try to retrieve data from 'x' and 'y' components
            # this only works if the data argument was used to fit the model
            if(!is.null(x <- object[["x"]]) && !is.null(y <- object$y)) {
                x <- removeIntercept(x)
                data <- data.frame(y, x)
            }
        }
        if(!haveDataArgument || is.null(data)) {
            # try to retrieve data from terms component
            data <- try(model.frame(object$terms), silent=TRUE)
            if(inherits(data, "try-error")) stop("model data not available")
        }
    }
    if(is.null(y <- object$y)) y <- model.response(data)
    ## call function cvFit() to perform cross-validation
    out <- cvFit(object, data=data, y=y, cost=cost, K=K, R=R, 
        foldType=foldType, folds=folds, costArgs=list(...), 
        envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


## MM and SDMD regression
#' @rdname repCV
#' @method repCV lmrob
#' @export
#' @import robustbase

repCV.lmrob <- function(object, cost=rtmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    # retrieve data from model fit
    if(is.null(data <- object$model)) {
        haveDataArgument <- !is.null(object$call$data)
        if(haveDataArgument) {
            # try to retrieve data from 'x' and 'y' components
            # this only works if the data argument was used to fit the model
            if(!is.null(x <- object[["x"]]) && !is.null(y <- object$y)) {
                x <- removeIntercept(x)
                data <- data.frame(y, x)
            }
        }
        if(!haveDataArgument || is.null(data)) {
            # try to retrieve data from terms component
            data <- try(model.frame(object$terms), silent=TRUE)
            if(inherits(data, "try-error")) stop("model data not available")
        }
    }
    if(is.null(y <- object$y)) y <- model.response(data)
    ## call function cvFit() to perform cross-validation
    out <- cvFit(object, data=data, y=y, cost=cost, K=K, R=R, 
        foldType=foldType, folds=folds, costArgs=list(...), 
        envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


## LTS regression
#' @rdname repCV
#' @method repCV lts
#' @export
#' @import robustbase

repCV.lts <- function(object, cost = rtmspe, K = 5, R = 1, 
        foldType = c("random", "consecutive", "interleaved"), 
        folds = NULL, fit = c("reweighted", "raw", "both"), 
        seed = NULL, ...) {
    ## initializations
    matchedCall <- match.call()
    object <- object
    if(is.null(x <- object$X) || is.null(y <- object$Y)) {
        if(is.null(data <- object$model)) {
            if(is.null(x)) x <- try(model.matrix(object$terms), silent=TRUE)
            if(is.null(y)) y <- try(model.response(object$terms), silent=TRUE)
            if(inherits(x, "try-error") || inherits(y, "try-error")) {
                stop("model data not available")
            }
        } else {
            x <- model.matrix(object$terms, data)
            y <- model.response(data)
        }
    }
    # predictor matrix is stored with column for intercept (if any)
    x <- removeIntercept(x)
    ## prepare cross-validation
    # extract function call for model fit
    call <- object$call
    call[[1]] <- as.name("ltsReg")
    # if the model was fitted with formula method, 'formula' and 'data' 
    # arguments are removed from call and 'x' and 'y' are used instead
    call$formula <- NULL
    call$data <- NULL
    call$intercept <- object$intercept
    ## call function cvFit() to perform cross-validation
    out <- cvFit(call, x=x, y=y, cost=cost, K=K, R=R, foldType=foldType, 
        folds=folds, predictArgs=list(fit=fit), costArgs=list(...), 
        envir=parent.frame(), seed=seed)
    out$call <- matchedCall
    out
}


# ----------------------

## LS regression
#' @rdname repCV
#' @export
cvLm <- repCV.lm

## MM and SDMD regression
#' @rdname repCV
#' @export
#' @import robustbase
cvLmrob <- repCV.lmrob

## LTS regression
#' @rdname repCV
#' @export
#' @import robustbase
cvLts <- repCV.lts

# ----------------------

#' @S3method predict lts
#' @import robustbase

# there is no predict() method for "lts" objects in package 'robustbase'
predict.lts <- function(object, newdata, 
        fit = c("reweighted", "raw", "both"), ...) {
    ## initializations
    fit <- match.arg(fit)
    coef <- switch(fit,
        reweighted=coef(object),
        raw=object$raw.coefficients,
        both=cbind(reweighted=coef(object), raw=object$raw.coefficients))
    terms <- delete.response(object$terms)  # extract terms for model matrix
    if(missing(newdata) || is.null(newdata)) {
        if(is.null(newdata <- object$X)) {
            if(is.null(data <- object$model)) {
                newdata <- try(model.matrix(terms), silent=TRUE)
                if(inherits(newdata, "try-error")) {
                    stop("model data not available")
                }
            } else newdata <- model.matrix(terms, data)
        }
    } else {
        # interpret vector as row
        if(is.null(dim(newdata))) newdata <- t(newdata)
        # check dimensions if model was not specified with a formula, 
        # otherwise use the terms object to extract model matrix
        if(is.null(terms)) {
            newdata <- as.matrix(newdata)
            if(object$intercept) {
                # if model has an intercept, add a column of ones to the new 
                # data matrix (unless it already contains intercept column)
                newdata <- addIntercept(newdata, check=TRUE)
            }
            # check dimensions of new data
            p <- if(is.null(dim(coef))) length(coef) else nrow(coef)
            if(ncol(newdata) != p) {
                stop(sprintf("new data must have %d columns", p))
            }
        } else newdata <- model.matrix(terms, as.data.frame(newdata))
    }
    ## compute predictions
    # ensure that a vector is returned if only one fit is requested
    drop(newdata %*% coef)
}

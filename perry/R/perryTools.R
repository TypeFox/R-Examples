# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## compute the predictions

# generic function to be extensible
perryPredictions <- function(call, data = NULL, x = NULL, y, splits, 
        predictFun = predict, predictArgs = list(), names = NULL, 
        envir = parent.frame(), cl = NULL) {
    UseMethod("perryPredictions", splits)
}

# default method for built-in procedures
perryPredictions.default <- function(call, data = NULL, x = NULL, y, 
        splits, predictFun = predict, predictArgs = list(), names = NULL, 
        envir = parent.frame(), cl = NULL) {
    # initializations
    useParallel <- !is.null(cl)
    R <- splits$R
    # define an expression that obtains predictions in one replication
    if(is.null(data)) {
        if(is.null(names)) names <- c("x", "y")
        if(inherits(splits, "cvFolds")) 
            fun <- if(useParallel && R == 1) pcvXY else cvXY
        else if(inherits(splits, "randomSplits")) fun <- rsXY
        else if(inherits(splits, "bootSamples")) fun <- bootXY
        else stop("invalid data splits")
    } else {
        if(is.null(names)) names <- "data"
        if(inherits(splits, "cvFolds")) 
            fun <- if(useParallel && R == 1) pcvData else cvData
        else if(inherits(splits, "randomSplits")) fun <- rsData
        else if(inherits(splits, "bootSamples")) fun <- bootData
        else stop("invalid data splits")
    }
    # obtain list of predictions for all replications
    if(useParallel) {
        if(inherits(splits, "cvFolds") && R == 1) {
            lapply(seq_len(R), function(r) {
                    s <- getIndices(splits, r)
                    fun(s, call=call, data=data, x=x, y=y, 
                        predictFun=predictFun, predictArgs=predictArgs, 
                        names=names, envir=envir, cl=cl)
                })
        } else {
            parLapply(cl, seq_len(R), function(r) {
                    s <- getIndices(splits, r)
                    fun(s, call=call, data=data, x=x, y=y, 
                        predictFun=predictFun, predictArgs=predictArgs, 
                        names=names, envir=envir)
                })
        }
    } else {
        lapply(seq_len(R), function(r) {
                s <- getIndices(splits, r)
                fun(s, call=call, data=data, x=x, y=y, predictFun=predictFun, 
                    predictArgs=predictArgs, names=names, envir=envir)
            })
    }
}

# one replication of cross validation for functions that take the predictors 
# and the response as separate arguments
cvXY <- function(folds, call, data, x, y, predictFun, predictArgs, names, envir) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- lapply(folds, rsXY, call=call, x=x, y=y, predictFun=predictFun, 
        predictArgs=predictArgs, names=names, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# one replication of cross validation for functions that that have one argument 
# for all the data
cvData <- function(folds, call, data, x, y, predictFun, predictArgs, names, envir) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- lapply(folds, rsData, call=call, data=data, predictFun=predictFun, 
        predictArgs=predictArgs, names=names, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# one replication of cross validation via parallel computing for functions that 
# take the predictors and the response as separate arguments
pcvXY <- function(folds, call, data, x, y, predictFun, predictArgs, names, envir, cl) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    # parLapply() already has an argument 'x', so don't use the argument names 
    # for the data or it will throw an error
    tmp <- parLapply(cl, folds, rsXY, call, data, x, y, 
        predictFun=predictFun, predictArgs=predictArgs, 
        names=names, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# one replication of cross validation via parallel computing for functions that 
# that have one argument for all the data
pcvData <- function(folds, call, data, x, y, predictFun, predictArgs, names, envir, cl) {
    # fit the model leaving each block out and obtain predictions for the 
    # left-out block
    tmp <- parLapply(cl, folds, rsData, call=call, data=data, 
        predictFun=predictFun, predictArgs=predictArgs, 
        names=names, envir=envir)
    # instead of collecting the results from the folds in the original order 
    # of the observations, they are simply stacked on top of each other
    combineData(tmp)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that take the predictors and the response as separate arguments
rsXY <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, -i)
    call[[names[2]]] <- dataSubset(y, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(x, i), args=predictArgs)
}

# fit the model for the training data and obtain predictions of the test data
# for functions that have one argument for all the data
rsData <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, -i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(data, i), args=predictArgs)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that take the predictors and the 
# response as separate arguments
bootXY <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names[1]]] <- dataSubset(x, i)
    call[[names[2]]] <- dataSubset(y, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(x, -i), args=predictArgs)
}

# fit the model for the bootstrap sample and obtain predictions for the 
# out-of-bag observations for functions that have one argument for all the data
bootData <- function(i, call, data, x, y, predictFun, predictArgs, names, envir) {
    # plug training data into function call
    call[[names]] <- dataSubset(data, i)
    # evaluate function call in supplied environment to make sure 
    # that other arguments are evaluated correctly
    fit <- eval(call, envir)
    # predict response for test data
    doCall(predictFun, fit, dataSubset(data, -i), args=predictArgs)
}


## estimate the prediction loss

# generic function to be extensible
perryCost <- function(splits, y, yHat, cost, costArgs = list()) {
    UseMethod("perryCost")
}

# default method for built-in procedures
perryCost.default <- function(splits, y, yHat, cost, costArgs = list()) {
    # initializations
    if(inherits(splits, "cvFolds")) {
        # function to compute prediction loss for all observations
        # response needs to be re-oredered according to cross-validation folds
        fun <- function(r, keepSE) {
            s <- unlist(getIndices(splits, r), use.names=FALSE)
            computeCost(cost, dataSubset(y, s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else if(inherits(splits, "randomSplits")) { 
        # function to compute prediction loss for test data
        fun <- function(r, keepSE) {
            s <- getIndices(splits, r)
            computeCost(cost, dataSubset(y, s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else if(inherits(splits, "bootSamples")) {
        # function to compute prediction loss for out-of-bag observations
        fun <- function(r, keepSE) {
            s <- getIndices(splits, r)
            computeCost(cost, dataSubset(y, -s), yHat[[r]], 
                args=costArgs, keepSE=keepSE)
        }
    } else stop("invalid data splits")
    # compute prediction loss
    R <- splits$R
    boot0.632 <- inherits(splits, "bootSamples") && splits$type == "0.632"
    if(R == 1) {
        # keep standard error if returned by prediction loss function, 
        # otherwise set it to NA
        pe <- fun(r=1, keepSE=!boot0.632)
        if(!is.list(pe)) {
            se <- rep.int(NA, length(pe))
            names(se) <- names(pe)
            pe <- list(pe=pe, se=se)
        }
        pe <- addNames(pe)
    } else {
        # discard standard error if returned by prediction loss function and 
        # recompute it from replications
        pe <- lapply(seq_len(R), fun, keepSE=FALSE)  # for each replication
        pe <- addNames(combineData(pe, drop=FALSE))  # combine results
    }
    # if requested, compute the 0.632 bootstrap estimator
    if(boot0.632) {
        # check if the fitted values from the model using all observations 
        # are available
        if(is.null(yHat <- splits$yHat)) 
            stop("fitted values to compute the apparent error are not available")
        else ae <- computeCost(cost, y, yHat, args=costArgs, keepSE=FALSE)
        # compute 0.632 estimator as combination of apparent error and 
        # out-of-bag prediction error
        if(R == 1) pe <- 0.632 * pe + 0.368 * ae
        else pe <- sweep(0.632 * pe, 2, 0.368 * ae, "+", check.margin=FALSE)
    }
    # aggregate results for more than one replication
    if(R > 1) {
        reps <- pe
        pe <- list(pe=apply(reps, 2, mean), se=apply(reps, 2, sd), reps=reps)
    }
    # return list
    pe
}

# compute cost function for predicted values from one replication
computeCost <- function(fun, y, yHat, args = list(), keepSE = TRUE) {
    # if the response is a vector and the predicted values are a matrix, 
    # compute the cost for each column of the matrix of predictions
    if(is.null(dim(y)) && !is.null(dim(yHat))) {
        pe <- apply(yHat, 2, function(x) doCall(fun, y, x, args=args))
        if(is.list(pe)) {
            # cost function returns list of prediction error and standard error
            if(keepSE) 
                pe <- list(pe=sapply(pe, "[[", 1), se=sapply(pe, "[[", 2))
            else pe <- sapply(pe, "[[", 1)
        }
    } else {
        pe <- doCall(fun, y, yHat, args=args)
        if(is.list(pe)) {
            pe <- if(keepSE) list(pe=pe[[1]], se=pe[[2]]) else pe[[1]]
        }
    }
    pe
}

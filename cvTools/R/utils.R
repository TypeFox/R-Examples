# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## utilities for cross-validation functions

# retrieve the number of observations
nobs.default <- function(object, ...) {
    n <- nrow(object)                   # matrix or data.frame
    if(is.null(n)) n <- length(object)  # vector
    n
}

# retrieve data subsets
dataSubset <- function(x, i, drop = FALSE) {
    if(is.null(dim(x))) {
        x[i]
    } else x[i, , drop=FALSE]
}

# replace data subsets
"dataSubset<-" <- function(x, i, value) {
    if(is.null(dim(x))) {
        x[i] <- value
    } else x[i, ] <- value
    x
}

# combine data (used for predictions from CV folds)
combineData <- function(x) {
    if(is.null(dim(x[[1]]))) {
        unlist(x)
    } else do.call("rbind", x)
}

## call a function by either
# 1) simply evaluating a supplied function for the basic arguments if there are
#    no additional arguments in list format
# 2) evaluating a supplied function with 'do.call' if there are additional 
#    arguments in list format
doCall <- function(fun, ..., args) {
    if(length(args) == 0) {
        fun(...)
    } else do.call(fun, c(list(...), args))
}

### call the supplied function for selecting the respective best model and 
### modify the multiplication factor of the standard error as necessary
#findBest <- function(fun, x, se, seFactor = NULL, returnSeFactor = TRUE) {
#    formalArgs <- formals(fun)
#    argNames <- names(formalArgs)
#    if(length(argNames) == 1) {
#        # only one argument for the prediction errors, call the supplied 
#        # function directly
#        best <- sapply(x[, -1, drop=FALSE], fun)
#        seFactor <- NA
#    } else {
#        # more than one argument: first should expect the prediction errors and 
#        # second the corresponding standard errors
#        if("seFactor" %in% argNames) {
#            # function contains argument for the multiplication factor of the 
#            # standard error
#            # if this is not is supplied, store the default value of the
#            # supplied function
#            # otherwise include the supplied value as extra argument in a list
#            # (this way it can easily be extended to take additional arguments)
#            if(is.null(seFactor)) {
#                seFactor <- formalArgs$seFactor
#                args <- list()
#            } else args <- list(seFactor=seFactor)
#        } else {
#            # function does not contain argument for the multiplication factor 
#            # of the standard error
#            args <- list()
#            seFactor <- NA
#        }
#        # call the supplied function via doCall()
#        best <- sapply(names(x)[-1], 
#            function(j) doCall(fun, x[, j], se[, j], args=args))
#    }
#    if(returnSeFactor) {
#        # return list with the first component containing the respective best 
#        # model and the second component containing the multiplication factor 
#        # of the standard error
#        list(best=best, seFactor=seFactor)
#    } else best
#}

# default names
defaultCvNames <- function(p) {
    if(p == 1) {
        "CV"
    } else if(p > 0) {
        paste("CV", seq_len(p), sep="")
    } else character()
}
defaultFitNames <- function(m) {
    if(m == 1) {
        "Fit"
    } else if(m > 0) {
        paste("Fit", seq_len(m), sep="")
    } else character()
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
    if(!check || all(is.na(match(c("Intercept","(Intercept)"), colnames(x))))) {
        cbind("(Intercept)"=rep.int(1, nrow(x)), x)
    } else x
}

## remove intercept column from design matrix
removeIntercept <- function(x, pos) {
    if(missing(pos)) {
        pos <- match(c("Intercept","(Intercept)"), colnames(x), nomatch = 0)
        if(any(pos > 0)) x[, -pos, drop=FALSE] else x
    } else x[, -pos, drop=FALSE]
}

# ----------------------

## utilities for plot functions

# get formula for plot functions
getFormula <- function(left, right, conditional = NULL) {
    if(is.null(conditional)) {
        as.formula(paste(left, "~", right))
    } else as.formula(paste(left, "~", right, "|", conditional))
}


# get data in the correct format for lattice graphics
getLatticeData <- function(x, ...) UseMethod("getLatticeData")

getLatticeData.cv <- function(x, select = NULL, reps = TRUE, 
        seFactor = NA, ...) {
    # extract subset of models
    x <- subset(x, select=select)
    if(reps) {
        CV <- x$reps
        if(is.null(CV)) {
            stop("replications not available")
        } else CV <- as.data.frame(CV)
    } else {
        CV <- as.data.frame(t(x$cv))
    }
    # stack selected results on top of each other
    cvName <- defaultCvNames(1)
    cvNames <- cvNames(x)
    ncv <- ncv(x)
    n <- nrow(CV)
    if(ncv == 0) {
        # return data frame of NAs if column is not selected
        CV[, cvName] <- rep.int(as.numeric(NA), n)
        if(!reps) SE <- as.numeric(NA)
    } else {
        if(!isTRUE(cvNames == cvName)) {
            CV <- lapply(cvNames, 
                function(j) data.frame(Name=rep.int(j, n), CV=CV[, j]))
            CV <- do.call(rbind, CV)
            names(CV) <- c("Name", cvName)
        }
        if(!reps) SE <- unname(x$se)
    }
    # return data
    if(reps) {
        CV
    } else {
        halflength <- seFactor * SE
        lower <- CV[, cvName] - halflength
        upper <- CV[, cvName] + halflength
        list(CV=CV, lower=lower, upper=upper)
    }
}

getLatticeData.cvSelect <- function(x, subset = NULL, select = NULL, 
        reps = TRUE, seFactor = x$seFactor, numericAsFactor = FALSE, ...) {
    # extract subset of models
    x <- subset(x, subset=subset, select=select)
    fits <- fits(x)
    if(reps) {
        CV <- x$reps
        if(is.null(CV)) stop("replications not available")
    } else {
        CV <- x$cv
        SE <- x$se
    }
    # ensure that models are shown in the correct order and drop unused levels
    # ensure that correct values are shown for a numeric tuning parameter
    if(numericAsFactor && is.double(CV$Fit)) {
        CV$Fit <- factor(shingle(CV$Fit), levels=fits)
        if(!reps) SE$Fit <- factor(shingle(SE$Fit), levels=fits)
    } else if(numericAsFactor || !is.numeric(CV$Fit)) {
        CV$Fit <- factor(CV$Fit, levels=fits)
        if(!reps) SE$Fit <- factor(SE$Fit, levels=fits)
    }
    # stack selected results on top of each other
    cvName <- defaultCvNames(1)
    cvNames <- cvNames(x)
    nfits <- nfits(x)
    ncv <- ncv(x)
    n <- nrow(CV)
    if(nfits == 0) {
        # no models selected: no column for grouping
        if(isTRUE(cvNames == cvName) || ncv == 0) {
            # return data frame without column for conditional plots and one NA 
            CV <- data.frame(as.numeric(NA))
            names(CV) <- cvName
            if(!reps) SE <- as.numeric(NA)
        } else {
            # return data frame with column for conditional plots and NA values
            CV <- data.frame(cvNames, rep.int(as.numeric(NA), ncv))
            names(CV) <- c("Name", cvName)
            if(!reps) SE <- rep.int(as.numeric(NA), ncv)
        }
    } else {
        # include column for grouping
        if(ncv == 0) {
            # no results selected: no column for conditional plots and NA values
            CV <- CV[, "Fit", drop=FALSE]
            CV[, cvName] <- rep.int(as.numeric(NA), n)
            if(!reps) SE <- rep.int(as.numeric(NA), nfits)
        } else {
            # no column for conditional plots if there is only one column of 
            # results with default name
            if(isTRUE(cvNames == cvName)) {
                if(!reps) SE <- SE[, cvName]
            } else {
                CVFit <- CV[, "Fit", drop=FALSE]
                CV <- lapply(cvNames, 
                    function(j) cbind(CVFit, Name=rep.int(j, n), CV=CV[, j]))
                CV <- do.call(rbind, CV)
                names(CV) <- c("Fit", "Name", cvName)
                if(!reps) SE <- unlist(SE[, cvNames], use.names=FALSE)
            }
        }
    }
    # return data
    if(reps) {
        CV
    } else {
        if(is.null(seFactor)) seFactor <- NA
        halflength <- seFactor * SE
        lower <- CV[, cvName] - halflength
        upper <- CV[, cvName] + halflength
        list(CV=CV, lower=lower, upper=upper)
    }
}

getLatticeData.cvTuning <- function(x, ...) {
    # adjust column specifying the model in case of only one tuning parameter
    if(ncol(x$tuning) == 1) fits(x) <- x$tuning[, 1]
    # call method for class "cvSelect"
    CV <- getLatticeData.cvSelect(x, ...)
    # return data
    CV
}

# Author: Harry Southworth
# Date: 2013-1-9
# Purpose: Split out worker functions to perform boring tasks for extreme
#          value modelling functions.
#
###########################################################################

texmexMethod <-
    # Take character string passed by user and coerce to standard format
function(method){
    method <- casefold(method)
    if (method %in% c("o", "opt", "optim", "optimize", "optimise")){
        method <- "o"
    }
    else if (method %in% c("s", "sim", "simulate")){
        method <- "s"
    }
    else if (method %in% c("b", "bs", "boot", "bootstrap")){
        method <- "b"
    }
    else {
        stop("method should be either 'optimize', 'simulate' or 'bootstrap'")
    }
    method
}

texmexPrior <-
    # Take input passed by user and coerce to standard format
function(prior, penalty, method, pp){
    prior <- casefold(prior)
    penalty <- casefold(penalty)

    # Deal with default situation and do MLE
    if (length(penalty) == 0 & is.null(pp) & method != 's'){
      prior <- "none"
    }
    else if (length(penalty) > 0){
      if (penalty != prior){
        prior <- penalty
      }
    }

    if (method == 's' & !is.element(prior, c('gaussian', 'cauchy'))){
      stop('Only gaussian or cauchy prior can be used when simulating from posterior.')
    }
    prior
}

texmexTrace <-
    # Get tracing frequency for optimizer and for Markov chain
function(trace, method){
    if (method == "o"){
      if (!is.null(trace)){ # trace provided by user
          otrace <- trace
      }
      else {
          otrace <- 0
      }
    } # Close if (method == "o"
    else{ # method == "s"
      otrace <- 0
      if (is.null(trace)){
         trace <- 10000
      }
    }
    c(otrace, trace)
}

texmexPrepareData <-
    # Get design matrices
function(y, data, params){
    D <- vector('list', length=length(params))
    if (!is.null(data)){
        y <- formula(paste(y, "~ 1"))
        y <- model.response(model.frame(y, data=data))

        for (i in 1:length(params)){
          D[[i]] <- model.matrix(params[[i]], data)
        }
    } # Close if(!is.null(data
    else {
        for (i in 1:length(params)){
            if (length(as.character(params[[i]])) == 2 & as.character(params[[i]])[2] == "1"){
                D[[i]] <- matrix(ncol = 1, rep(1, length(y)))
            }
            else {
                D[[i]] <- model.matrix(params[[i]])
            }
        } # Close for
    } # Close else

    names(D) <- names(params)

    list(y=y, D=D)
}

texmexThresholdData <- function(threshold, data){
    # Need to subset design matrices on y > th, so do those
    # first, then threshold y

    for (i in 1:length(data$D)){
        data$D[[i]] <- data$D[[i]][data$y > threshold, , drop=FALSE]
    }

    data$y <- data$y[data$y > threshold]
    if (length(data$y) == 0){
      stop("No observations above the threshold.")
    }

    data
}

texmexPriorParameters <-
    # Pre-process prior distribution parameters
function(prior, priorParameters, data){

    # Get total number of parameters
    nc <- sum(sapply(data$D, ncol))

    if (prior %in% c("quadratic", "gaussian")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, nc), diag(rep(10^4, nc)))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Gaussian prior or quadratic penalty, priorParameters should be a list of length 2, the second element of which should be a symmetric (covariance) matrix")
        }
    }
    else if (prior %in% c("lasso", "l1", "laplace")) {
        if (is.null(priorParameters)) {
            priorParameters <- list(rep(0, nc), diag(rep(10^(-4), nc)))
        }
        if (length(priorParameters) != 2 | !is.list(priorParameters)) {
            stop("For Laplace prior or L1 or Lasso penalty, priorParameters should be a list of length 2, the second element of which should be a diagonal (precision) matrix")
        }
        if (!is.matrix(priorParameters[[2]])) {
            priorParameters[[2]] <- diag(rep(priorParameters[[2]], nc))
        }
        if (!all(priorParameters[[2]] == diag(diag(priorParameters[[2]])))) {
            warning("some off-diagonal elements of the covariance are non-zero. Only the diagonal is used in penalization")
        }
    }

    #### If priorParameters given but of wrong dimension, kill
    if (!is.null(priorParameters)) {
        if (length(priorParameters[[1]]) != nc) {
            stop("wrong number of parameters in prior (doesn't match parameter formulae)")
        }
        else if (length(diag(priorParameters[[2]])) != nc) {
            stop("wrong dimension of prior covariance (doesn't match parameter formulae)")
        }
    }
    priorParameters
}

findFormulae <-
    # Find formulae in a call
function(call,...){
    wh <- sapply(call, function(x){ try(class(eval(x)), silent=TRUE) })
    wh <- names(wh)[wh == 'formula']
    if (length(wh) > 0){
        res <- as.list(call[wh])
        for(i in 1:length(res)) res[[i]] <- eval(res[[i]]) # don't use sapply - can't handle "..." arguments
    }
    else { res <- NULL }
    res
}

texmexParameters <- function(call, fam, ...){
    # Create intercept formula for every parameter
    mp <- lapply(fam$param, function(x) ~1)
    names(mp) <- fam$param

    # Splice in parameters from function call
    p <- findFormulae(call, ...)
    mp[names(p)] <- p

    mp
}

texmexGetParam <- function(data, co){
    if (is.vector(co)){ co <- matrix(co, nrow=1) }

    mend <- cumsum(unlist(lapply(data, ncol)))
    mstart <- c(1, mend+1)[-(length(mend) + 1)]

    param <- lapply(1:length(mend), function(i, m, start, end){
                                        m[,start[i]:end[i],drop=FALSE]
                                    },
                    m=co, start=mstart, end=mend)
    param
}

texmexPst <- function(msg="",Family){
  paste(msg,", Family = ",Family$name)
}

texmexGetXlevels <- function(fo, data){
  # Get all variable names used on RHSs of formulae
  getVars <- function(fo) { all.vars(update(fo, 0~.)) }
  allVars <- unique(unlist(lapply(fo, getVars)))

  # Get rid of variables not in data, get classes, then get rid of non-factors
  data <- data[, allVars, drop=FALSE]
  classes <- sapply(data, class)
  data <- data[, classes %in% c("factor", "ordered", "character"), drop=FALSE]
  
  data[classes == "character"] <- lapply(data[classes == "character"], as.factor)
  
  # Get a single named list containing all levels
  xlevels <- lapply(data, levels)
  
  # Split it by formula
  res <- lapply(fo, getVars)
  lapply(res, function(X, wh) wh[names(wh) %in% X], wh=xlevels)
}

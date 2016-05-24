# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Resampling-based prediction error for tuning parameter selection
#' 
#' Select tuning parameters of a model by estimating the respective prediction 
#' errors via (repeated) \eqn{K}-fold cross-validation, (repeated) random 
#' splitting (also known as random subsampling or Monte Carlo 
#' cross-validation), or the bootstrap.  It is thereby possible to supply a 
#' model fitting function or an unevaluated function call to a model fitting 
#' function.
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
#' @aliases coef.perryTuning fitted.perryTuning predict.perryTuning 
#' print.perryTuning residuals.perryTuning
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
#' vector of values can be supplied.  The prediction error is then estimated 
#' for all possible combinations of tuning parameter values.
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
#' @param final  a logical indicating whether to fit the final model with the 
#' optimal combination of tuning parameters.
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
#' @return 
#' If \code{tuning} is an empty list, \code{\link{perryFit}} is called to 
#' return an object of class \code{"perry"}.
#' 
#' Otherwise an object of class \code{"perryTuning"} (which inherits from class 
#' \code{"perrySelect"}) with the following components is returned:
#' @returnItem pe  a data frame containing the estimated prediction errors for 
#' all combinations of tuning parameter values.  In case of more than one 
#' replication, those are average values over all replications.
#' @returnItem se  a data frame containing the estimated standard errors of the 
#' prediction loss for all combinations of tuning parameter values.
#' @returnItem reps  a data frame containing the estimated prediction errors 
#' from all replications for all combinations of tuning parameter values.  This 
#' is only returned in case of more than one replication.
#' @returnItem splits  an object giving the data splits used to estimate the 
#' prediction error.
#' @returnItem y  the response.
#' @returnItem yHat  a list containing the predicted values for all 
#' combinations of tuning parameter values.  Each list component is again a 
#' list containing the corresponding predicted values from all replications.
#' @returnItem best  an integer vector giving the indices of the optimal 
#' combinations of tuning parameters.
#' @returnItem selectBest  a character string specifying the criterion used for 
#' selecting the best model.
#' @returnItem seFactor  a numeric value giving the multiplication factor of 
#' the standard error used for the selection of the best model.
#' @returnItem tuning  a data frame containing the grid of tuning parameter 
#' values for which the prediction error was estimated.
#' @returnItem finalModel  the final model fit with the optimal combination of 
#' tuning parameters.  This is only returned if argument \code{final} is 
#' \code{TRUE}.
#' @returnItem call  the matched function call.
#' 
#' @note 
#' The same data splits are used for all combinations of tuning parameter 
#' values for maximum comparability.
#' 
#' If a final model with the optimal combination of tuning parameters is 
#' computed, class \code{"perryTuning"} inherits the \code{coef()}, 
#' \code{fitted()}, \code{predict()} and \code{residuals()} methods from 
#' its component \code{finalModel}.
#' 
#' @author Andreas Alfons
#' 
#' @references 
#' Hastie, T., Tibshirani, R. and Friedman, J. (2009) \emph{The Elements of 
#' Statistical Learning: Data Mining, Inference, and Prediction}.  Springer, 
#' 2nd edition.
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{cvFolds}}, \code{\link{randomSplits}}, 
#' \code{\link{bootSamples}}, \code{\link{cost}}
#' 
#' @example inst/doc/examples/example-perryTuning.R
#' 
#' @keywords utilities
#' 
#' @export

perryTuning <- function(object, ...) UseMethod("perryTuning")


#' @rdname perryTuning
#' @method perryTuning function
#' @export

perryTuning.function <- function(object, formula, data = NULL, x = NULL, y, 
                                 tuning = list(), args = list(), 
                                 splits = foldControl(), predictFun = predict, 
                                 predictArgs = list(), cost = rmspe, 
                                 costArgs = list(), 
                                 selectBest = c("min", "hastie"), seFactor = 1, 
                                 final = FALSE, names = NULL, 
                                 envir = parent.frame(), ncores = 1, cl = NULL, 
                                 seed = NULL, ...) {
  ## initializations
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("perryTuning")
  # set up unevaluated function call
  final <- isTRUE(final)
  call <- as.call(c(if(final) substitute(object) else object, args))
  # check formula and data
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
  out <- perryTuning(call, data=data, x=x, y=y, tuning=tuning, splits=splits, 
                     predictFun=predictFun, predictArgs=predictArgs, cost=cost, 
                     costArgs=costArgs, selectBest=selectBest, 
                     seFactor=seFactor, final=final, names=names, envir=envir, 
                     ncores=ncores, cl=cl, seed=seed, ...)
  out$call <- matchedCall
  out
}


#' @rdname perryTuning
#' @method perryTuning call
#' @export

perryTuning.call <- function(object, data = NULL, x = NULL, y, tuning = list(), 
                             splits = foldControl(), predictFun = predict, 
                             predictArgs = list(), cost = rmspe, 
                             costArgs = list(), selectBest = c("min", "hastie"), 
                             seFactor = 1, final = FALSE, names = NULL, 
                             envir = parent.frame(), ncores = 1, cl = NULL, 
                             seed = NULL, ...) {
  ## initializations
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("perryTuning")
  n <- nobs(y)
  if(is.null(data)) {
    sx <- "x"
    nx <- nobs(x)
  } else {
    sx <- "data"
    nx <- nobs(data)
  }
  if(!isTRUE(n == nx)) stop(sprintf("'%s' must have %d observations", sx, nx))
  # create all combinations of tuning parameters
  tuning <- do.call(expand.grid, tuning)
  nTuning <- nrow(tuning)
  pTuning <- ncol(tuning)
  if(nTuning == 0 || pTuning == 0) {
    # use function perryFit() if no tuning parameters are supplied
    out <- perryFit(object, data, x, y, splits=splits, predictFun=predictFun, 
                    predictArgs=predictArgs, cost=cost, costArgs=costArgs, 
                    names=names, envir=envir, ncores=ncores, cl=cl, seed=seed)
    return(out)
  }
  if(!is.null(seed)) set.seed(seed)
  ## compute data splits
  if(hasMethod("perrySplits", class(splits))) splits <- perrySplits(n, splits)
  # set up parallel computing if requested
  haveCl <- inherits(cl, "cluster")
  if(haveCl) haveNcores <- FALSE
  else {
    if(is.na(ncores)) ncores <- detectCores()  # use all available cores
    if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
      ncores <- 1  # use default value
      warning("invalid value of 'ncores'; using default value")
    } else ncores <- as.integer(ncores)
    if(nTuning == 1) {
      R <- splits$R
      if(inherits(splits, "cvFolds") && R == 1) {
        ncores <- min(ncores, splits$K)
      } else ncores <- min(ncores, R)
    } else ncores <- min(ncores, nTuning)
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
  ## compute the predictions for each combination of tuning parameters
  fits <- seq_len(nTuning)
  tn <- names(tuning)
  if(useParallel) {
    # set seed of the random number stream
    if(!is.null(seed)) clusterSetRNGStream(cl, iseed=seed)
    else if(haveNcores) clusterSetRNGStream(cl)
    if(nTuning == 1) {
      yHat <- lapply(fits, function(i) {
        # add tuning parameters to function call
        for(j in seq_len(pTuning)) object[[tn[j]]] <- tuning[i, j]
        # compute the predictions
        perryPredictions(object, data, x, y, splits=splits, 
                         predictFun=predictFun, predictArgs=predictArgs, 
                         names=names, envir=envir, cl=cl)
      })
    } else {
      yHat <- parLapply(cl, fits, function(i) {
        # add tuning parameters to function call
        for(j in seq_len(pTuning)) object[[tn[j]]] <- tuning[i, j]
        # compute the predictions
        perryPredictions(object, data, x, y, splits=splits, 
                         predictFun=predictFun, predictArgs=predictArgs, 
                         names=names, envir=envir)
      })
    }
  } else {
    yHat <- lapply(fits, function(i) {
      # add tuning parameters to function call
      for(j in seq_len(pTuning)) object[[tn[j]]] <- tuning[i, j]
      # compute the predictions
      perryPredictions(object, data, x, y, splits=splits, 
                       predictFun=predictFun, predictArgs=predictArgs, 
                       names=names, envir=envir)
    })
  }
  ## estimate the prediction loss for each combination of tuning parameters
  pe <- lapply(yHat, function(yHat) {
    perryCost(splits, y, yHat, cost=cost, costArgs=costArgs)
  })
  pe <- combineResults(pe, fits=fits)
  ## select optimal tuning parameters
  best <- selectBest(pe$pe, pe$se, method=selectBest, seFactor=seFactor)
  ## compute final model if requested
  final <- isTRUE(final)
  if(final) {
    # plug data into function call
    if(is.null(data)) {
      if(is.null(names)) names <- c("x", "y")
      object[[names[1]]] <- substitute(x)
      object[[names[2]]] <- substitute(y)
    } else {
      if(is.null(names)) names <- "data"
      object[[names]] <- substitute(data)
    }
    # add optimal combination of tuning parameters to function call
    i <- best$best
    for(j in seq_len(pTuning)) object[[tn[j]]] <- unique(tuning[i, j])
    # evaluate function call to compute final model
    # use this environment since data are added with substitute()
    finalModel <- try(eval(object))
    if(inherits(finalModel, "try-error")) {
      final <- FALSE
      warn <- gsub("Error in", "In", finalModel)
      warning(warn, call.=FALSE)
    }
  }
  ## construct return object
  names(yHat) <- fits
  pe <- c(pe, list(splits=splits, y=y, yHat=yHat), best, list(tuning=tuning))
  if(final) pe$finalModel <- finalModel
  pe$call <- matchedCall
  class(pe) <- c("perryTuning", "perrySelect")
  pe
}

# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' (Robust) groupwise least angle regression
#' 
#' (Robustly) sequence groups of candidate predictors according to their 
#' predictive content and find the optimal model along the sequence.
#' 
#' @aliases print.grplars
#' 
#' @param formula  a formula describing the full model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{grplars} or \code{rgrplars} is called.
#' @param x  a matrix or data frame containing the candidate predictors.
#' @param y  a numeric vector containing the response.
#' @param sMax  an integer giving the number of predictor groups to be 
#' sequenced.  If it is \code{NA} (the default), predictor groups are sequenced 
#' as long as there are twice as many observations as expected predictor 
#' variables (number of predictor groups times the average number of predictor 
#' variables per group).
#' @param assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @param centerFun  a function to compute a robust estimate for the center 
#' (defaults to \code{\link[stats]{median}}).
#' @param scaleFun  a function to compute a robust estimate for the scale 
#' (defaults to \code{\link[stats]{mad}}).
#' @param regFun  a function to compute robust linear regressions that can be 
#' interpreted as weighted least squares (defaults to 
#' \code{\link[robustbase]{lmrob}}).
#' @param regArgs  a list of arguments to be passed to \code{regFun}.
#' @param combine  a character string specifying how to combine the data 
#' cleaning weights from the robust regressions with each predictor group.  
#' Possible values are \code{"min"} for taking the minimum weight for each 
#' observation, \code{"euclidean"} for weights based on Euclidean distances 
#' of the multivariate set of standardized residuals (i.e., multivariate 
#' winsorization of the standardized residuals assuming independence), or 
#' \code{"mahalanobis"} for weights based on Mahalanobis distances of the 
#' multivariate set of standardized residuals (i.e., multivariate winsorization 
#' of the standardized residuals).
#' @param const  numeric; tuning constant for multivariate winsorization to be 
#' used in the initial corralation estimates based on adjusted univariate 
#' winsorization (defaults to 2).
#' @param prob  numeric; probability for the quantile of the 
#' \eqn{\chi^{2}}{chi-squared} distribution to be used in multivariate 
#' winsorization (defaults to 0.95).
#' @param fit  a logical indicating whether to fit submodels along the sequence 
#' (\code{TRUE}, the default) or to simply return the sequence (\code{FALSE}).
#' @param s  an integer vector of length two giving the first and last 
#' step along the sequence for which to compute submodels.  The default 
#' is to start with a model containing only an intercept (step 0) and 
#' iteratively add all groups along the sequence (step \code{sMax}).  If 
#' the second element is \code{NA}, predictor groups are added to the 
#' model as long as there are twice as many observations as predictor 
#' variables.  If only one value is supplied, it is recycled.
#' @param crit  a character string specifying the optimality criterion to be 
#' used for selecting the final model.  Possible values are \code{"BIC"} for 
#' the Bayes information criterion and \code{"PE"} for resampling-based 
#' prediction error estimation.
#' @param splits  an object giving data splits to be used for prediction error 
#' estimation (see \code{\link[perry]{perry}}).
#' @param cost  a cost function measuring prediction loss (see 
#' \code{\link[perry]{perry}} for some requirements).  The 
#' default is to use the root trimmed mean squared prediction error for a 
#' robust fit and the root mean squared prediction error otherwise (see 
#' \code{\link[perry]{cost}}).
#' @param costArgs  a list of additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' @param selectBest,seFactor  arguments specifying a criterion for selecting 
#' the best model (see \code{\link[perry]{perrySelect}}).  The default is to 
#' use a one-standard-error rule.
#' @param ncores  a positive integer giving the number of processor cores to be 
#' used for parallel computing (the default is 1 for no parallelization).  If 
#' this is set to \code{NA}, all available processor cores are used.  For 
#' obtaining the data cleaning weights, for fitting models along the sequence 
#' and for prediction error estimation, parallel computing is implemented on 
#' the \R level using package \pkg{parallel}.  Otherwise parallel computing for 
#' some of of the more computer-intensive computations in the sequencing step 
#' is implemented on the C++ level via OpenMP (\url{http://openmp.org/}).
#' @param cl  a \pkg{parallel} cluster for parallel computing as generated by 
#' \code{\link[parallel]{makeCluster}}.  This is preferred over \code{ncores} 
#' for tasks that are parallelized on the \R level, in which case \code{ncores} 
#' is only used for tasks that are parallelized on the C++ level.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).  This is useful because many robust regression 
#' functions (including \code{\link[robustbase]{lmrob}}) involve randomness, 
#' or for prediction error estimation.  On parallel \R worker processes, random 
#' number streams are used and the seed is set via 
#' \code{\link{clusterSetRNGStream}}.
#' @param model  a logical indicating whether the model data should be included 
#' in the returned object.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return 
#' If \code{fit} is \code{FALSE}, an integer vector containing the indices of 
#' the sequenced predictor groups.
#'  
#' Else if \code{crit} is \code{"PE"}, an object of class 
#' \code{"perrySeqModel"} (inheriting from classes \code{"perryTuning"}, 
#' see \code{\link[perry]{perryTuning}}).  It contains information on the 
#' prediction error criterion, and includes the final model as component 
#' \code{finalModel}.
#' 
#' Otherwise an object of class \code{"grplars"} (inheriting from class 
#' \code{"seqModel"}) with the following components:
#' @returnItem active  an integer vector containing the sequence of predictor 
#' groups.
#' @returnItem s  an integer vector containing the steps for which submodels 
#' along the sequence have been computed.
#' @returnItem coefficients  a numeric matrix in which each column contains the 
#' regression coefficients of the corresponding submodel along the sequence.
#' @returnItem fitted.values  a numeric matrix in which each column contains 
#' the fitted values of the corresponding submodel along the sequence.
#' @returnItem residuals  a numeric matrix in which each column contains 
#' the residuals of the corresponding submodel along the sequence.
#' @returnItem df  an integer vector containing the degrees of freedom of the 
#' submodels along the sequence (i.e., the number of estimated coefficients).
#' @returnItem robust  a logical indicating whether a robust fit was computed.
#' @returnItem scale  a numeric vector giving the robust residual scale 
#' estimates for the submodels along the sequence (only returned for a robust 
#' fit).
#' @returnItem crit  an object of class \code{"bicSelect"} containing the BIC 
#' values and indicating the final model (only returned if argument \code{crit} 
#' is \code{"BIC"} and argument \code{s} indicates more than one step along the 
#' sequence).
#' @returnItem muX  a numeric vector containing the center estimates of the 
#' predictor variables.
#' @returnItem sigmaX  a numeric vector containing the scale estimates of the 
#' predictor variables.
#' @returnItem muY  numeric; the center estimate of the response.
#' @returnItem sigmaY  numeric; the scale estimate of the response.
#' @returnItem x  the matrix of candidate predictors (if \code{model} is 
#' \code{TRUE}).
#' @returnItem y  the response (if \code{model} is \code{TRUE}).
#' @returnItem assign  an integer vector giving the predictor group to which 
#' each predictor variable belongs.
#' @returnItem w  a numeric vector giving the data cleaning weights (only 
#' returned for a robust fit).
#' @returnItem call  the matched function call.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link[=coef.seqModel]{coef}}, 
#' \code{\link[=fitted.seqModel]{fitted}}, 
#' \code{\link[=plot.seqModel]{plot}}, 
#' \code{\link[=predict.seqModel]{predict}}, 
#' \code{\link[=residuals.seqModel]{residuals}}, 
#' \code{\link[robustbase]{lmrob}}
#' 
#' @example inst/doc/examples/example-rgrplars.R
#' 
#' @keywords regression robust
#' 
#' @export

grplars <- function(x, ...) UseMethod("grplars")


#' @rdname grplars
#' @method grplars formula
#' @export

grplars.formula <- function(formula, data, ...) {
  ## initializations
  call <- match.call()  # get function call
  call[[1]] <- as.name("grplars")
  # prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 1  # ensure model with intercept
  if(is.empty.model(mt)) stop("empty model")
  # extract response and candidate predictors from model frame
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  ## call default method
  out <- grplars.default(x, y, ...)
  if(inherits(out, "grplars")) {
    out$call <- call  # add call to return object
    out$terms <- mt   # add model terms to return object
  }
  out
}


#' @rdname grplars
#' @method grplars data.frame
#' @export

grplars.data.frame <- function(x, y, ...) {
  ## initializations
  call <- match.call()  # get function call
  call[[1]] <- as.name("grplars")
  x <- model.matrix(~ ., data=x)   # convert data.frame to design matrix
  ## call default method
  out <- grplars.default(x, y, ...)
  if(inherits(out, "grplars")) out$call <- call  # add call to return object
  out
}


#' @rdname grplars
#' @method grplars default
#' @export

grplars.default <- function(x, y, sMax = NA, assign, fit = TRUE, s = c(0, sMax), 
                            crit = c("BIC", "PE"), splits = foldControl(), 
                            cost = rmspe, costArgs = list(), 
                            selectBest = c("hastie", "min"), seFactor = 1, 
                            ncores = 1, cl = NULL, seed = NULL, model = TRUE, 
                            ...) {
  ## initializations
  call <- match.call()  # get function call
  call[[1]] <- as.name("grplars")
  # if argument 'assign' is not supplied, check if predictor matrix has 
  # attribute "assign" (as generated by model.matrix())
  if(missing(assign)) {
    assign <- attr(x, "assign")
    if(is.null(assign)) {
      # take each predictor to be its own group
      assign <- seq_len(ncol(x))
    }
  }
  # check if the predictor matrix contains column for intercept and 
  # remove it if necessary
  if(isTRUE(assign[1] == 0)) {
    x <- removeIntercept(x)
    assign <- assign[-1]
  }
  ## call fit function with classical functions for center, scale, 
  ## correlation and regression
  grplarsFit(x, y, sMax=sMax, assign=assign, robust=FALSE, centerFun=mean, 
             scaleFun=sd, fit=fit, s=s, crit=crit, splits=splits, cost=cost, 
             costArgs=costArgs, selectBest=selectBest, seFactor=seFactor, 
             ncores=ncores, cl=cl, seed=seed, model=model, call=call)
}


#' @rdname grplars
#' @export

rgrplars <- function(x, ...) UseMethod("rgrplars")


#' @rdname grplars
#' @method rgrplars formula
#' @export

rgrplars.formula <- function(formula, data, ...) {
  ## initializations
  call <- match.call()  # get function call
  call[[1]] <- as.name("rgrplars")
  # prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  attr(mt, "intercept") <- 1  # ensure model with intercept
  if(is.empty.model(mt)) stop("empty model")
  # extract response and candidate predictors from model frame
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  ## call default method
  out <- rgrplars.default(x, y, ...)
  if(inherits(out, "grplars")) {
    out$call <- call  # add call to return object
    out$terms <- mt   # add model terms to return object
  }
  out
}


#' @rdname grplars
#' @method rgrplars data.frame
#' @export

rgrplars.data.frame <- function(x, y, ...) {
  ## initializations
  call <- match.call()  # get function call
  call[[1]] <- as.name("rgrplars")
  x <- model.matrix(~ ., data=x)   # convert data.frame to design matrix
  ## call default method
  out <- rgrplars.default(x, y, ...)
  if(inherits(out, "grplars")) out$call <- call  # add call to return object
  out
}


#' @rdname grplars
#' @method rgrplars default
#' @export

rgrplars.default <- function(x, y, sMax = NA, assign, centerFun = median, 
                             scaleFun = mad, regFun = lmrob, regArgs = list(), 
                             combine = c("min", "euclidean", "mahalanobis"), 
                             const = 2, prob = 0.95, fit = TRUE, 
                             s = c(0, sMax), crit = c("BIC", "PE"), 
                             splits = foldControl(), cost = rtmspe, 
                             costArgs = list(), selectBest = c("hastie", "min"), 
                             seFactor = 1, ncores = 1, cl = NULL, seed = NULL, 
                             model = TRUE, ...) {
  ## initializations
  call <- match.call()  # get function call
  call[[1]] <- as.name("rgrplars")
  # if argument 'assign' is not supplied, check if predictor matrix has 
  # attribute "assign" (as generated by model.matrix())
  if(missing(assign)) {
    assign <- attr(x, "assign")
    if(is.null(assign)) {
      # take each predictor to be its own group
      assign <- seq_len(ncol(x))
    }
  }
  # check if the predictor matrix contains column for intercept and 
  # remove it if necessary
  if(isTRUE(assign[1] == 0)) {
    x <- removeIntercept(x)
    assign <- assign[-1]
  }
  ## call fit function with supplied functions for center, scale, 
  ## correlation and regression
  grplarsFit(x, y, sMax=sMax, assign=assign, robust=TRUE, centerFun=centerFun, 
             scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, 
             combine=combine, const=const, prob=prob, fit=fit, s=s, crit=crit, 
             splits=splits, cost=cost, costArgs=costArgs, selectBest=selectBest, 
             seFactor=seFactor, ncores=ncores, cl=cl, seed=seed, model=model, 
             call=call)
}


## fit function that allows to specify functions for center, scale, correlation 
## and regression
grplarsFit <- function(x, y, sMax = NA, assign, robust = FALSE, 
                       centerFun = mean, scaleFun = sd, 
                       regFun = lm.fit, regArgs = list(), 
                       combine = c("min", "euclidean", "mahalanobis"), 
                       const = 2, prob = 0.95, fit = TRUE, s = c(0, sMax), 
                       crit = c("BIC", "PE"), splits = foldControl(), 
                       cost = rmspe, costArgs = list(), 
                       selectBest = c("hastie", "min"), seFactor = 1, 
                       ncores = 1, cl = NULL, seed = NULL, model = TRUE, 
                       call = NULL) {
  ## initializations
  n <- length(y)
  x <- addColnames(as.matrix(x))
  if(nrow(x) != n) stop(sprintf("'x' must have %d rows", n))
  ## call workhorse function
  grouplars(x, y, sMax=sMax, assign=assign, robust=robust, centerFun=centerFun, 
            scaleFun=scaleFun, regFun=regFun, regArgs=regArgs, combine=combine, 
            const=const, prob=prob, fit=fit, s=s, crit=crit, splits=splits, 
            cost=cost, costArgs=costArgs, selectBest=selectBest, 
            seFactor=seFactor, ncores=ncores, cl=cl, seed=seed, model=model, 
            call=call)
}

# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

#' Sparse least trimmed squares regression
#' 
#' Compute least trimmed squares regression with an \eqn{L_{1}}{L1} penalty on 
#' the regression coefficients, which allows for sparse model estimates.
#' 
#' @aliases print.sparseLTS
#' 
#' @param formula  a formula describing the model.
#' @param data  an optional data frame, list or environment (or object coercible 
#' to a data frame by \code{\link{as.data.frame}}) containing the variables in 
#' the model.  If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which 
#' \code{sparseLTS} is called.
#' @param x  a numeric matrix containing the predictor variables.
#' @param y  a numeric vector containing the response variable.
#' @param lambda  a numeric vector of non-negative values to be used as penalty 
#' parameter.
#' @param mode  a character string specifying the type of penalty parameter.  If 
#' \code{"lambda"}, \code{lambda} gives the grid of values for the penalty 
#' parameter directly.  If \code{"fraction"}, the smallest value of the penalty 
#' parameter that sets all coefficients to 0 is first estimated based on 
#' bivariate winsorization, then \code{lambda} gives the fractions of that 
#' estimate to be used (hence all values of \code{lambda} should be in the 
#' interval [0,1] in that case).
#' @param alpha  a numeric value giving the percentage of the residuals for 
#' which the \eqn{L_{1}}{L1} penalized sum of squares should be minimized (the 
#' default is 0.75).
#' @param normalize  a logical indicating whether the predictor variables 
#' should be normalized to have unit \eqn{L_{2}}{L2} norm (the default is 
#' \code{TRUE}).  Note that normalization is performed on the subsamples 
#' rather than the full data set.
#' @param intercept  a logical indicating whether a constant term should be 
#' included in the model (the default is \code{TRUE}).
#' @param nsamp  a numeric vector giving the number of subsamples to be used in 
#' the two phases of the algorithm.  The first element gives the number of 
#' initial subsamples to be used.  The second element gives the number of 
#' subsamples to keep after the first phase of \code{ncstep} C-steps.  For 
#' those remaining subsets, additional C-steps are performed until 
#' convergence.  The default is to first perform \code{ncstep} C-steps on 500 
#' initial subsamples, and then to keep the 10 subsamples with the lowest value 
#' of the objective function for additional C-steps until convergence.
#' @param initial  a character string specifying the type of initial subsamples 
#' to be used.  If \code{"sparse"}, the lasso fit given by three randomly 
#' selected data points is first computed.  The corresponding initial subsample 
#' is then formed by the fraction \code{alpha} of data points with the smallest 
#' squared residuals.  Note that this is optimal from a robustness point of 
#' view, as the probability of including an outlier in the initial lasso fit is 
#' minimized.  If \code{"hyperplane"}, a hyperplane through \eqn{p} randomly 
#' selected data points is first computed, where \eqn{p} denotes the number of 
#' variables.  The corresponding initial subsample is then again formed by the 
#' fraction \code{alpha} of data points with the smallest squared residuals.  
#' Note that this cannot be applied if \eqn{p} is larger than the number of 
#' observations.  Nevertheless, the probability of including an outlier 
#' increases with increasing dimension \eqn{p}.  If \code{"random"}, the 
#' initial subsamples are given by a fraction \code{alpha} of randomly 
#' selected data points.  Note that this leads to the largest probability of 
#' including an outlier.
#' @param ncstep  a positive integer giving the number of C-steps to perform on 
#' all subsamples in the first phase of the algorithm (the default is to 
#' perform two C-steps).
#' @param use.correction  currently ignored.  Small sample correction factors 
#' may be added in the future.
#' @param tol  a small positive numeric value giving the tolerance for 
#' convergence.
#' @param eps  a small positive numeric value used to determine whether the 
#' variability within a variable is too small (an effective zero).
#' @param use.Gram  a logical indicating whether the Gram matrix of the 
#' explanatory variables should be precomputed in the lasso fits on the 
#' subsamples.  If the number of variables is large, computation may be faster 
#' when this is set to \code{FALSE}.  The default is to use \code{TRUE} if the 
#' number of variables is smaller than the number of observations in the 
#' subsamples and smaller than 100, and \code{FALSE} otherwise.
#' @param crit  a character string specifying the optimality criterion to be 
#' used for selecting the final model.  Possible values are \code{"BIC"} for 
#' the Bayes information criterion and \code{"PE"} for resampling-based 
#' prediction error estimation.
#' @param splits  an object giving data splits to be used for prediction error 
#' estimation (see \code{\link[perry]{perryTuning}}).
#' @param cost  a cost function measuring prediction loss (see 
#' \code{\link[perry]{perryTuning}} for some requirements).  The 
#' default is to use the root trimmed mean squared prediction error 
#' (see \code{\link[perry]{cost}}).
#' @param costArgs  a list of additional arguments to be passed to the 
#' prediction loss function \code{cost}.
#' @param selectBest,seFactor  arguments specifying a criterion for selecting 
#' the best model (see \code{\link[perry]{perryTuning}}).  The default is to 
#' use a one-standard-error rule.
#' @param ncores  a positive integer giving the number of processor cores to be 
#' used for parallel computing (the default is 1 for no parallelization).  If 
#' this is set to \code{NA}, all available processor cores are used.  For 
#' prediction error estimation, parallel computing is implemented on the \R 
#' level using package \pkg{parallel}.  Otherwise parallel computing is 
#' implemented on the C++ level  via OpenMP (\url{http://openmp.org/}).
#' @param cl  a \pkg{parallel} cluster for parallel computing as generated by 
#' \code{\link[parallel]{makeCluster}}.  This is preferred over \code{ncores} 
#' for prediction error estimation, in which case \code{ncores} is only used on 
#' the C++ level for computing the final model.
#' @param seed  optional initial seed for the random number generator (see 
#' \code{\link{.Random.seed}}).  On parallel \R worker processes for prediction 
#' error estimation, random number streams are used and the seed is set via 
#' \code{\link{clusterSetRNGStream}}.
#' @param model  a logical indicating whether the data \code{x} and \code{y} 
#' should be added to the return object.  If \code{intercept} is \code{TRUE}, 
#' a column of ones is added to \code{x} to account for the intercept.
#' @param \dots  additional arguments to be passed down.
#' 
#' @return 
#' If \code{crit} is \code{"PE"}, an object of class \code{"perrySparseLTS"} 
#' (inheriting from class \code{"perryTuning"}, see 
#' \code{\link[perry]{perryTuning}}).  It contains information on the 
#' prediction error criterion, and includes the final model with the optimal 
#' tuning paramter as component \code{finalModel}.
#' 
#' Otherwise an object of class \code{"sparseLTS"} with the following 
#' components:
#' @returnItem lambda  a numeric vector giving the values of the penalty 
#' parameter.
#' @returnItem best  an integer vector or matrix containing the respective best 
#' subsets of \eqn{h} observations found and used for computing the raw 
#' estimates.
#' @returnItem objective  a numeric vector giving the respective values of the 
#' sparse LTS objective function, i.e., the \eqn{L_{1}}{L1} penalized sums of 
#' the \eqn{h} smallest squared residuals from the raw fits.
#' @returnItem coefficients  a numeric vector or matrix containing the 
#' respective coefficient estimates from the reweighted fits.
#' @returnItem fitted.values  a numeric vector or matrix containing the 
#' respective fitted values of the response from the reweighted fits.
#' @returnItem residuals  a numeric vector or matrix containing the 
#' respective residuals from the reweighted fits.
#' @returnItem center  a numeric vector giving the robust center estimates of 
#' the corresponding reweighted residuals.
#' @returnItem scale  a numeric vector giving the robust scale estimates of the 
#' corresponding reweighted residuals.
#' @returnItem cnp2  a numeric vector giving the respective consistency factors 
#' applied to the scale estimates of the reweighted residuals.
#' @returnItem wt  an integer vector or matrix containing binary weights that 
#' indicate outliers from the respective reweighted fits, i.e., the weights are 
#' \eqn{1} for observations with reasonably small reweighted residuals and 
#' \eqn{0} for observations with large reweighted residuals.
#' @returnItem df  an integer vector giving the respective degrees of freedom 
#' of the obtained reweighted model fits, i.e., the number of nonzero 
#' coefficient estimates.
#' @returnItem intercept  a logical indicating whether the model includes a 
#' constant term.
#' @returnItem alpha  a numeric value giving the percentage of the residuals for 
#' which the \eqn{L_{1}}{L1} penalized sum of squares was minimized.
#' @returnItem quan  the number \eqn{h} of observations used to compute the raw 
#' estimates.
#' @returnItem raw.coefficients  a numeric vector or matrix containing the 
#' respective coefficient estimates from the raw fits.  
#' @returnItem raw.fitted.values  a numeric vector or matrix containing the 
#' respective fitted values of the response from the raw fits.
#' @returnItem raw.residuals  a numeric vector or matrix containing the 
#' respective residuals from the raw fits.
#' @returnItem raw.center  a numeric vector giving the robust center estimates 
#' of the corresponding raw residuals.
#' @returnItem raw.scale  a numeric vector giving the robust scale estimates of 
#' the corresponding raw residuals.
#' @returnItem raw.cnp2  a numeric value giving the consistency factor applied 
#' to the scale estimate of the raw residuals.
#' @returnItem raw.wt  an integer vector or matrix containing binary weights 
#' that indicate outliers from the respective raw fits, i.e., the weights used 
#' for the reweighted fits.
#' @returnItem crit  an object of class \code{"bicSelect"} containing the BIC 
#' values and indicating the final model (only returned if argument \code{crit} 
#' is \code{"BIC"} and argument \code{lambda} contains more than one value for 
#' the penalty parameter).
#' @returnItem x  the predictor matrix (if \code{model} is \code{TRUE}).
#' @returnItem y  the response variable (if \code{model} is \code{TRUE}).
#' @returnItem call  the matched function call.
#' 
#' @note Package \pkg{robustHD} has a built-in back end for sparse least 
#' trimmed squares using the C++ library Armadillo.  Another back end is 
#' available through package \pkg{sparseLTSEigen}, which uses the C++ library 
#' Eigen.  The latter is faster, currently does not work on 32-bit \R for 
#' Windows.
#' 
#' For both C++ back ends, parallel computing is implemented via OpenMP 
#' (\url{http://openmp.org/}).
#' 
#' @author Andreas Alfons
#' 
#' @references
#' Alfons, A., Croux, C. and Gelper, S. (2013) Sparse least trimmed squares 
#' regression for analyzing high-dimensional large data sets. \emph{The Annals 
#' of Applied Statistics}, \bold{7}(1), 226--248.
#' 
#' @seealso \code{\link[=coef.sparseLTS]{coef}}, 
#' \code{\link[=fitted.sparseLTS]{fitted}}, 
#' \code{\link[=plot.sparseLTS]{plot}}, 
#' \code{\link[=predict.sparseLTS]{predict}}, 
#' \code{\link[=residuals.sparseLTS]{residuals}}, 
#' \code{\link{wt}}, \code{\link[robustbase]{ltsReg}}
#' 
#' @example inst/doc/examples/example-sparseLTS.R
#' 
#' @keywords regression robust
#' 
#' @export 
#' @import parallel
#' @import perry

sparseLTS <- function(x, ...) UseMethod("sparseLTS")


#' @rdname sparseLTS
#' @method sparseLTS formula
#' @export

sparseLTS.formula <- function(formula, data, ...) {
  ## get function call
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("sparseLTS")
  ## prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if(is.empty.model(mt)) stop("empty model")
  ## extract response and candidate predictors from model frame
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  ## call default method
  # check if the specified model contains an intercept
  # if so, remove the column for intercept and use 'intercept=TRUE'
  # otherwise use 'intercept=FALSE'
  whichIntercept <- match("(Intercept)", colnames(x), nomatch = 0)
  intercept <- whichIntercept > 0
  if(intercept) {
    x <- x[, -whichIntercept, drop = FALSE]
    fit <- sparseLTS(x, y, intercept=TRUE, ...)
  } else fit <- sparseLTS(x, y, intercept=FALSE, ...)
  # return results
  fit$call <- matchedCall  # add call to return object
  fit$terms <- mt          # add model terms to return object
  fit
}


#' @rdname sparseLTS
#' @method sparseLTS default
#' @export

sparseLTS.default <- function(x, y, lambda, mode = c("lambda", "fraction"), 
                              alpha = 0.75, normalize = TRUE, intercept = TRUE, 
                              nsamp = c(500, 10), 
                              initial = c("sparse", "hyperplane", "random"), 
                              ncstep = 2, use.correction = TRUE, 
                              tol = .Machine$double.eps^0.5, 
                              eps = .Machine$double.eps, use.Gram, 
                              crit = c("BIC", "PE"), splits = foldControl(), 
                              cost = rtmspe, costArgs = list(), 
                              selectBest = c("hastie", "min"), 
                              seFactor = 1, ncores = 1, cl = NULL, 
                              seed = NULL, model = TRUE, ...) {
  ## initializations
  matchedCall <- match.call()
  matchedCall[[1]] <- as.name("sparseLTS")
  n <- length(y)
  x <- addColnames(as.matrix(x))
  d <- dim(x)
  if(!isTRUE(n == d[1])) stop(sprintf("'x' must have %d rows", n))
  if(missing(lambda)) {
    # if penalty parameter is not supplied, use a small fraction of a 
    # robust estimate of the smallest value that sets all coefficients 
    # to zero
    lambda <- 0.05
    mode <- "fraction"
  } else {
    # otherwise check the supplied penalty parameter
    if(!is.numeric(lambda) || length(lambda) == 0 || any(!is.finite(lambda))) {
      stop("missing or invalid value of 'lambda'")
    }
    if(any(negative <- lambda < 0)) {
      lambda[negative] <- 0
      warning("negative value for 'lambda', using no penalization")
    }
    lambda <- sort.int(unique(lambda), decreasing=TRUE)
    mode <- match.arg(mode)
  }
  if(length(lambda) == 1) crit <- "none" 
  else crit <- match.arg(crit)
  normalize <- isTRUE(normalize)
  intercept <- isTRUE(intercept)
  if(mode == "fraction" && any(lambda > 0)) { 
    # fraction of a robust estimate of the smallest value for the penalty 
    # parameter that sets all coefficients to zero (based on bivariate 
    # winsorization)
    if(crit == "PE") frac <- lambda
    lambda <- lambda * lambda0(x, y, normalize=normalize, intercept=intercept, 
                               tol=tol, eps=eps, ...)
  }
  alpha <- rep(alpha, length.out=1)
  if(!isTRUE(is.numeric(alpha) && 0.5 <= alpha && alpha <= 1)) {
    stop("'alpha' must be between 0.5 and 1")
  }
  nsamp <- rep(nsamp, length.out=2)
  if(!is.numeric(nsamp) || any(!is.finite(nsamp))) {
    nsamp <- formals$nsamp()
    warning("missing or infinite values in 'nsamp'; using default values")
  } else nsamp <- as.integer(nsamp)
  ncstep <- rep(as.integer(ncstep), length.out=1)
  if(!is.numeric(ncstep) || !is.finite(ncstep)) {
    ncstep <- formals()$ncstep
    warning("missing or infinite value of 'ncstep'; using default value")
  } else ncstep <- as.integer(ncstep)
  tol <- rep(tol, length.out=1)
  if(!is.numeric(tol) || !is.finite(tol)) {
    tol <- formals()$tol
    warning("missing or infinite value of 'tol'; using default value")
  }
  eps <- rep(eps, length.out=1)
  if(!is.numeric(eps) || !is.finite(eps)) {
    eps <- formals()$eps
    warning("missing or infinite value of 'eps'; using default value")
  }
  h <- ceiling(alpha*n)  # subset sizes are different from 'ltsReg'
  if(missing(use.Gram)) use.Gram <- d[2] >= min(h, 100)
  else use.Gram <- isTRUE(use.Gram)
  ncores <- rep(ncores, length.out=1)
  if(is.na(ncores)) ncores <- detectCores()  # use all available cores
  if(!is.numeric(ncores) || is.infinite(ncores) || ncores < 1) {
    ncores <- 1  # use default value
    warning("invalid value of 'ncores'; using default value")
  } else ncores <- as.integer(ncores)
  
  ## compute sparse LTS
  if(crit == "PE") {
    # set up function call to be passed to perryTuning()
    remove <- c("x", "y", "lambda", "crit", "splits", "cost", "costArgs", 
                "selectBest", "seFactor", "ncores", "cl", "seed")
    remove <- match(remove, names(matchedCall), nomatch=0)
    call <- matchedCall[-remove]
    # call function perryTuning() to perform prediction error estimation
    tuning <- list(lambda=if(mode == "fraction") frac else lambda)
    selectBest <- match.arg(selectBest)
    fit <- perryTuning(call, x=x, y=y, tuning=tuning, splits=splits, 
                       predictArgs=list(fit="both"), cost=cost, 
                       costArgs=costArgs, selectBest=selectBest, 
                       seFactor=seFactor, ncores=ncores, cl=cl, 
                       seed=seed)
    # fit final model
    lambdaOpt <- unique(lambda[fit$best])
    call$x <- matchedCall$x
    call$y <- matchedCall$y
    call$lambda <- lambdaOpt
    call$mode <- NULL
    call$ncores <- matchedCall$ncores
    finalModel <- eval(call)
    # if optimal tuning parameter is different for reweighted and raw fit, 
    # add information that indicates optimal values
    if(length(lambdaOpt)) {
      crit <- list(best=c(reweighted=1, raw=2))
      class(crit) <- "fitSelect"
      finalModel$crit <- crit
    }
    # modify results
    if(mode == "fraction" && any(lambda > 0)) fit$tuning$lambda <- lambda
    fit$finalModel <- finalModel
    class(fit) <- c("perrySparseLTS", class(fit))
  } else {
    if(h < n) {
      initial <- match.arg(initial)
      if(h < d[2] && initial == "hyperplane") initial <- "sparse"
      if(!is.null(seed)) set.seed(seed)  # set seed of random number generator
      
      ## the same initial subsets are used for all values of lambda, except when 
      ## the initial subsets are based on lasso solutions with 3 observations 
      ## (since in the latter case the initial subsets depend on lambda)
      subsets <- switch(initial, random=randomSubsets(n, h, nsamp[1]), 
                        hyperplane=hyperplaneSubsets(x, y, h, nsamp[1]))
      
      ## call internal function to obtain raw fits
      if(length(lambda) == 1) {
        fit <- fastSparseLTS(x=x, y=y, lambda=lambda, h=h, nsamp=nsamp, 
                             initial=subsets, normalize=normalize, 
                             intercept=intercept, ncstep=ncstep, tol=tol, 
                             eps=eps, use.Gram=use.Gram, ncores=ncores)
        fit$best <- sort.int(fit$best)
      } else {
        names(lambda) <- seq_along(lambda)
        fit <- lapply(lambda, fastSparseLTS, x=x, y=y, h=h, nsamp=nsamp, 
                      initial=subsets, normalize=normalize, intercept=intercept, 
                      ncstep=ncstep, tol=tol, eps=eps, use.Gram=use.Gram, 
                      ncores=ncores, drop=FALSE)
        names(fit) <- names(lambda)
        fit <- list(best=sapply(fit, function(x) sort.int(x$best)), 
                    coefficients=sapply(fit, "[[", "coefficients"), 
                    residuals=sapply(fit, "[[", "residuals"), 
                    objective=sapply(fit, "[[", "objective"), 
                    center=sapply(fit, "[[", "center"), 
                    scale=sapply(fit, "[[", "scale"))
      }
      ## compute consistency factor to correct scale estimate
      qn <- qnorm((h+n)/ (2*n))                         # required quantile
      cdelta <- 1 / sqrt(1 - (2*n)/(h/qn) * dnorm(qn))  # consistency factor
    } else {
      ## no trimming, compute lasso for raw fits
      if(length(lambda) == 1) {
        fit <- fastLasso(x, y, lambda=lambda, normalize=normalize, 
                         intercept=intercept, eps=eps, use.Gram=use.Gram, 
                         raw=TRUE)
      } else {
        names(lambda) <- seq_along(lambda)
        fit <- lapply(lambda, fastLasso, x=x, y=y, normalize=normalize, 
                      intercept=intercept, eps=eps, use.Gram=use.Gram, 
                      drop=FALSE, raw=TRUE)
        names(fit) <- names(lambda)
        fit <- list(best=sapply(fit, "[[", "best"), 
                    coefficients=sapply(fit, "[[", "coefficients"), 
                    residuals=sapply(fit, "[[", "residuals"), 
                    objective=sapply(fit, "[[", "objective"), 
                    center=sapply(fit, "[[", "center"), 
                    scale=sapply(fit, "[[", "scale"))
      }
      ## consistency factor is not necessary
      cdelta <- 1
    }
    
    ## find good observations
    q <- qnorm(0.9875)       # quantile of the normal distribution
    s <- fit$scale * cdelta  # corrected scale estimate
    if(length(lambda) == 1) ok <- abs((fit$residuals - fit$center)/s) <= q
    else ok <- abs(scale(fit$residuals, fit$center, s)) <= q
    ## convert to 0/1 weights identifying outliers
    raw.wt <- ok
    storage.mode(raw.wt) <- "integer"
    
    ## keep information on raw estimator
    raw.fit <- fit
    raw.cdelta <- cdelta
    raw.s <- s
    nOk <- if(length(lambda) == 1) sum(raw.wt) else colSums(raw.wt)
    
    ## compute reweighted estimator
    # compute lasso fits with good observations
    if(length(lambda) == 1) {
      fit <- fastLasso(x, y, lambda=lambda, subset=which(ok), 
                       normalize=normalize, intercept=intercept, 
                       eps=eps, use.Gram=use.Gram)
    } else {
      fit <- lapply(seq_along(lambda), function(i) {
        fastLasso(x, y, lambda=lambda[i], subset=which(ok[, i]), 
                  normalize=normalize, intercept=intercept, 
                  eps=eps, use.Gram=use.Gram, drop=FALSE)
      })
      names(fit) <- names(lambda)
      fit <- list(coefficients=sapply(fit, "[[", "coefficients"), 
                  fitted.values=sapply(fit, "[[", "fitted.values"), 
                  residuals=sapply(fit, "[[", "residuals"))
    }
    ## compute consistency factor
    qn <- qnorm((nOk+n)/ (2*n))  # quantile for consistency factor
    cdelta <- 1 / sqrt(1-(2*n)/(nOk/qn)*dnorm(qn))
    cdelta[nOk == n] <- 1
    ## compute 0/1 weights identifying outliers
    if(length(lambda) == 1) {
      # compute residual center and scale estimates
      center <- sum(raw.wt*fit$residuals)/nOk
      centeredResiduals <- fit$residuals - center
      s <- sqrt(sum(raw.wt*centeredResiduals^2)/(nOk-1)) * cdelta
      # compute outlier weights
      wt <- as.integer(abs(centeredResiduals/s) <= q)
    } else {
      # compute residual center and scale estimates
      center <- colSums(raw.wt*fit$residuals)/nOk
      centeredResiduals <- sweep(fit$residuals, 2, center, check.margin=FALSE)
      s <- sqrt(colSums(raw.wt*centeredResiduals^2)/(nOk-1)) * cdelta
      # compute outlier weights
      wt <- abs(sweep(centeredResiduals, 2, s, "/", check.margin=FALSE)) <= q
      storage.mode(wt) <- "integer"
    }
    
    ## construct return object
    if(intercept) x <- addIntercept(x)
    if(length(lambda) == 1) {
      df <- modelDf(fit$coefficients, tol)
      raw.df <- modelDf(raw.fit$coefficients, tol)
    } else {
      df <- apply(fit$coefficients, 2, modelDf, tol)
      raw.df <- apply(raw.fit$coefficients, 2, modelDf, tol)
    }
    fit <- list(lambda=lambda, best=raw.fit$best, objective=raw.fit$objective, 
                coefficients=copyNames(from=x, to=fit$coefficients), 
                fitted.values=copyNames(from=y, to=fit$fitted.values), 
                residuals=copyNames(from=y, to=fit$residuals), center=center, 
                scale=s, cnp2=cdelta, wt=copyNames(from=y, to=wt), df=df, 
                normalize=normalize, intercept=intercept, alpha=alpha, quan=h, 
                raw.coefficients=copyNames(from=x, to=raw.fit$coefficients), 
                raw.fitted.values=copyNames(from=y, to=y-raw.fit$residuals),
                raw.residuals=copyNames(from=y, to=raw.fit$residuals), 
                raw.center=raw.fit$center, raw.scale=raw.s, raw.cnp2=raw.cdelta,
                raw.wt=copyNames(from=y, to=raw.wt), raw.df=raw.df)
    class(fit) <- "sparseLTS"
    
    ## add information on the optimal model
    if(crit == "BIC") fit$crit <- bicSelect(fit, fit="both")
    
    ## add model data if requested
    if(isTRUE(model)) fit[c("x", "y")] <- list(x=x, y=y)
  }
  ## return results
  fit$call <- matchedCall
  fit
}


## internal function to compute raw sparse LTS
fastSparseLTS <- function(lambda, x, y, h, nsamp = c(500, 10), 
                          initial = NULL, normalize = TRUE, intercept = TRUE, 
                          ncstep = 2, tol = .Machine$double.eps^0.5, 
                          eps = .Machine$double.eps, use.Gram = TRUE, 
                          ncores = 1, drop = TRUE) {
  # check whether initial subsets based on lasso solutions need to be computed
  if(is.null(initial)) {
    initial <- sparseSubsets(x, y, lambda=lambda, h=h, nsamp=nsamp[1], 
                             normalize=normalize, intercept=intercept, 
                             eps=eps, use.Gram=use.Gram)
  }
  # call C++ function
  fit <- callBackend("R_fastSparseLTS", R_x=x, R_y=y, R_lambda=lambda, 
                     R_initial=initial, R_normalize=normalize, 
                     R_intercept=intercept, R_ncstep=ncstep, R_nkeep=nsamp[2], 
                     R_tol=tol, R_eps=eps, R_useGram=use.Gram, R_ncores=ncores)
  if(drop) {
    # drop the dimension of selected components
    which <- c("best", "coefficients", "residuals")
    fit[which] <- lapply(fit[which], drop)
  }
  fit
}

###########################################################################################

### internal functions taken from gamboostLSS 1.2-0
## helper functions
check <- function(what, what_char, names) {
  
  errormsg <- paste0(sQuote(what_char), " can be either a scalar, a (named) vector or a (named) list",
                     " of ", what_char, " values with same names as ",  sQuote("families"), "in ",
                     sQuote("boost_control"))
  
  if (is.list(what)) {
    if (is.null(names(what)) && length(what) == length(names))
      names(what) <- names
    if (!all(names(what) %in% names) ||
        length(unique(names(what))) != length(names))
      stop(errormsg)
    what <- what[names] ## sort in order of families
    what <- unlist(what)
  } else {
    if(length(what) != 1 && length(what) != length(names))
      stop(errormsg)
    if (length(what) == 1) {
      what <- rep(what, length(names))
      names(what) <- names
    } else {
      if (is.null(names(what)))
        names(what) <- names
      if (!all(names(what) %in% names))
        stop(errormsg)
      what <- what[names] ## sort in order of families
    }
  }
  
  return(what)
}


## helper function in a modified version based on mboost_2.2-3
## print trace of boosting iterations
do_trace <- function(current, mstart, risk,
                     linebreak = options("width")$width / 2, mstop = 1000) {
  current <- current - mstart
  if (current != mstop) {
    if ((current - 1) %/% linebreak == (current - 1) / linebreak) {
      mchr <- formatC(current + mstart, format = "d",
                      width = nchar(mstop) + 1, big.mark = "'")
      cat(paste("[", mchr, "] ",sep = ""))
    } else {
      if ((current %/% linebreak != current / linebreak)) {
        cat(".")
      } else {
        cat(" -- risk:", risk[current + mstart], "\n")
      }
    }
  } else {
    cat("\nFinal risk:", risk[current + mstart], "\n")
  }
}


## helper function from mboost 2.6-0
## mboost_intern(..., fun = "mboost_intern") gives error with one check
rescale_weights <- function(w) {
  if (max(abs(w - floor(w))) < sqrt(.Machine$double.eps))
    return(w)
  return(w / sum(w) * sum(w > 0))
}


## helper function copied from mboost_2.2-3
### check measurement scale of response for some losses
check_y_family <- function(y, family){
  ## <SB> check the response is in long format, otherwise is.numeric()
  ## as e.g. in StudentTMu()@check_y() gives an error
  family@check_y(as.vector(y)) 
  y
}



#################################################################################
#' Model-based Gradient Boosting for Functional GAMLSS
#' 
#' Function for fitting GAMLSS (generalized additive models for location, scale and shape) 
#' with functional data using component-wise gradient boosting, for details see 
#' Brockhaus et al. (2015a). 
#' 
#' @param formula a symbolic description of the model to be fit. 
#' If \code{formula} is a single formula, the same formula is used for all distribution parameters. 
#' \code{formula} can also be a (named) list, where each list element corresponds to one distribution 
#' parameter of the GAMLSS distribution. The names must be the same as in the \code{families}. 
#' @param timeformula one-sided formula for the expansion over the index of the response. 
#' For a functional response \eqn{Y_i(t)} typically \code{~bbs(t)} to obtain a smooth 
#' expansion of the effects along \code{t}. In the limiting case that \eqn{Y_i} is a scalar response
#' use \code{~bols(1)}, which sets up a base-learner for the scalar 1. 
#' Or you can use \code{timeformula=NULL}, then the scalar response is treated as scalar. 
#' Analogously to \code{formula}, \code{timeformula} can either be a one-sided formula or 
#' a named list of one-sided formulas. 
#' @param data a data frame or list containing the variables in the model.
#' @param families an object of class \code{families}. It can be either one of the pre-defined distributions 
#' that come along with the package \code{gamboostLSS} or a new distribution specified by the user 
#' (see \code{\link[gamboostLSS]{Families}} for details). 
#' Per default, the two-parametric \code{\link[gamboostLSS]{GaussianLSS}} family is used.
#' @param control  a list of parameters controlling the algorithm. 
#' For more details see \code{\link[mboost]{boost_control}}.  
#' @param weights does not work!
#' @param ... additional arguments passed to \code{\link[FDboost]{FDboost}}, 
#' including, \code{family} and \code{control}.
#' 
#' @details For details on the theory of GAMLSS see Rigby and Stasinopoulos (2005). 
#' \code{FDboostLSS} uses \code{FDboost} to fit the distibution parameters of a GAMLSS - 
#' a functional boosting model is fitted for each parameter. 
#' See \code{\link{FDboost}} for details on boosting functional regression models 
#' as introduced by Brockhaus et al. (2015b). 
#' See \code{\link[gamboostLSS]{mboostLSS}} for 
#' details on boosting of GAMLSS for scalar variables as introduced by Mayr et al. (2012). 
#' 
#' @return an object of class \code{FDboostLSS} that inherits from \code{mboostLSS}. 
#' The \code{FDboostLSS}-object is a named list containing one list entry per distribution parameter
#' and some attributes. The list is named like the parameters, e.g. mu and sigma, 
#' if the parameters mu and sigma are modelled. Each list-element is an object of class \code{FDboost}.  
#' 
#' @author Sarah Brockhaus 
#' 
#' @seealso Note that \code{FDboostLSS} calls \code{\link{FDboost}} directly.  
#' 
#' @keywords models, nonlinear 
#' 
#' @references 
#' Brockhaus, S. and Fuest, A. and Mayr, A. and Greven, S. (2015a): 
#' Functional regression models for location, scale and shape applied to stock returns. 
#' In: Friedl H. and Wagner H. (eds), 
#' Proceedings of the 30th International Workshop on Statistical Modelling: 117-122.  
#'
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015b). 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300.
#' 
#' Mayr, A., Fenske, N., Hofner, B., Kneib, T. and Schmid, M. (2012): 
#' Generalized additive models for location, scale and shape for high-dimensional 
#' data - a flexible approach based on boosting. 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 61(3), 403-427. 
#' 
#' Rigby, R. A. and D. M. Stasinopoulos (2005):  
#' Generalized additive models for location, scale and shape (with discussion). 
#' Journal of the Royal Statistical Society: Series C (Applied Statistics), 54(3), 507-554. 
#' 
#' @examples 
#' ########### simulate Gaussian scalar-on-function data
#' n <- 500 ## number of observations
#' G <- 120 ## number of observations per functional covariate
#' set.seed(123) ## ensure reproducibility
#' z <- runif(n) ## scalar covariate
#' z <- z - mean(z)
#' s <- seq(0, 1, l=G) ## index of functional covariate
#' ## generate functional covariate
#' if(require(splines)){
#'    x <- t(replicate(n, drop(bs(s, df = 5, int = TRUE) %*% runif(5, min = -1, max = 1))))
#' }else{
#'   x <- matrix(rnorm(n*G), ncol = G, nrow = n)
#' }
#' x <- scale(x, center = TRUE, scale = FALSE) ## center x per observation point
#' 
#' mu <- 2 + 0.5*z + (1/G*x) %*% sin(s*pi)*5 ## true functions for expectation
#' sigma <- exp(0.5*z - (1/G*x) %*% cos(s*pi)*2) ## for standard deviation
#' 
#' y <- rnorm(mean = mu, sd = sigma, n = n) ## draw respone y_i ~ N(mu_i, sigma_i)
#' 
#' ## save data as list containing s as well 
#' dat_list <- list(y = y, z = z, x = I(x), s = s)
#' 
#' ## model fit assuming Gaussian location scale model 
#' m_boost <- FDboostLSS(list(mu = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16), 
#'                            sigma = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16)), 
#'                            timeformula = NULL, data = dat_list)
#' summary(m_boost)
#' 
#' \dontrun{
#'  if(require(gamboostLSS)){
#'   ## find optimal number of boosting iterations on a grid in [1, 500]
#'   ## using 5-fold bootstrap
#'   grid <-  make.grid(c(mu = 500, sigma = 500), length.out = 10)
#'   ## takes some time, easy to parallelize on Linux
#'   set.seed(123) 
#'   cvr <- cvrisk(m_boost, folds = cv(model.weights(m_boost[[1]]), B = 5),
#'                 grid = grid, trace = FALSE)
#'   ## use model at optimal stopping iterations 
#'   m_boost <- m_boost[mstop(cvr)] ## [c(172, 63)]
#'    
#'   ## plot smooth effects of functional covariates
#'   par(mfrow = c(1,2))
#'   plot(m_boost$mu, which = 2, ylim = c(0,5))
#'   lines(s, sin(s*pi)*5, col = 3, lwd = 2)
#'   plot(m_boost$sigma, which = 2, ylim = c(-2.5,2.5))
#'   lines(s, -cos(s*pi)*2, col = 3, lwd = 2)
#'  }
#' }
#' @export
## function that calls FDboost for each distribution parameter
FDboostLSS <- function(formula, timeformula, data = list(), families = GaussianLSS(),
                       control = boost_control(), weights = NULL, ...){
  
  cl <- match.call()
  if(is.null(cl$families)) cl$families <- families
  
  ## warnings for functional response are irrelevant for scalar response  
  if( !is.null(timeformula) && timeformula != ~bols(1) ){
    message("No smooth offsets over time are used, just global scalar offsets.")
    message("No integration weights are used to compute the loss for the functional response.")
  }
  
  fit <- FDboostLSS_fit(formula = formula, timeformula = timeformula, 
                        data = data, families = families,
                        control = control, weights = weights, ...,
                        fun = FDboost, funchar = "FDboost", call = cl)
  return(fit)
}

###########################################################################################

### work horse for fitting FDboostLSS models
### based on code of mboostLSS_fit() from package gamboostLSS 1.2-0
FDboostLSS_fit <- function(formula, timeformula, data = list(), families = GaussianLSS(),
                           control = boost_control(), weights = NULL,
                           fun = FDboost, funchar = "FDboost", call = NULL, ...){
  
  if (length(families) == 0)
    stop(sQuote("families"), " not specified")
  
  if ("offset" %in% names(list(...)))
    stop("Do not use argument ", sQuote("offset"),
         ". Please specify offsets via families")
  ### Use mu in "families" to specify offset in mu-family etc.
  
  if (is.list(formula)){
    if (!all(names(formula) %in% names(families)) ||
        length(unique(names(formula))) != length(names(families)))
      stop(sQuote("formula"), " can be either a formula or a named list",
           " of formulas with same names as ",  sQuote("families"), ".")
    ynames <- sapply(formula, function(fm) as.character(fm[[2]]))
    if (length(unique(ynames)) > 1)
      warning("responses differ for the components")
    response <- lapply(formula, function(fm)
      eval(as.expression(fm[[2]]), envir = data,
           enclos = environment(fm)))
    #response <- eval(as.expression(formula[[1]][[2]]), envir = data,
    #                 enclos = environment(formula[[1]]))
  } else {
    response <- eval(as.expression(formula[[2]]), envir = data,
                     enclos = environment(formula))
    tmp <- vector("list", length = length(families))
    names(tmp) <- names(families)
    for (i in 1:length(tmp))
      tmp[[i]] <- formula
    formula <- tmp
  }
  
  if (is.list(timeformula)){
    if (!all(names(timeformula) %in% names(families)) ||
        length(unique(names(timeformula))) != length(names(families)))
      stop(sQuote("timeformula"), " can be either a one-sided formula or a named list",
           " of timeformulas with same names as ",  sQuote("families"), ".")
  } else {
    tmp <- vector("list", length = length(families))
    names(tmp) <- names(families)
    for (i in 1:length(tmp))
      tmp[i] <- list(timeformula)
    timeformula <- tmp
  }
  mstop <- mstoparg <- control$mstop
  control$mstop <- 1
  mstop <- check(mstop, "mstop", names(families))
  
  nu <- control$nu
  nu <- check(nu, "nu", names(families))
  
  if (is.list(control$risk) || is.list(control$center) || is.list(control$trace))
    stop(sQuote("risk"),", ", sQuote("center"), " and ", sQuote("trace") ,
         " cannot be lists in ", sQuote("boost_control"))
  
  trace <- control$trace
  control$trace <- FALSE
  
  ### <FIXME> SB: is the dealing with the weights correct??
  w <- weights 
  if (is.null(weights)){
    if (!is.list(response)) {
      weights <- rep.int(1, NROW(response))
      
      ### <SB> expand weights if the response is a matrix (functional response)
      if( !is.null(dim(response)) && !any(dim(response)==1) ){
        weights <- rep.int(weights, ncol(response))
      }
      
    } else {
      weights <- rep.int(1, NROW(response[[1]]))
      
      ### <SB> expand weights if the response is a matrix (functional response)
      if( !is.null(dim(response[[1]])) && !any(dim(response[[1]])==1) ){
        weights <- rep.int(weights, ncol(response[[1]]))
      }
      
    }
  } ## weights=rep(1, ncol(response[[j]])*nrow(response[[j]]) )
  
  ## weights <- mboost_intern(weights, fun = "rescale_weights")
  weights <- rescale_weights(weights)
  
  fit <- vector("list", length = length(families))
  names(fit) <- names(families)
  
  mods <- 1:length(fit)
  
  offset <- vector("list", length = length(mods))
  names(offset) <- names(families)
  for (j in mods){
    if (!is.list(response)) {
      response <- check_y_family(response, families[[j]])
      offset[[j]] <- families[[j]]@offset(y = c(response), w = weights)
    } else {
      response[[j]] <- check_y_family(response[[j]], families[[j]])
      offset[[j]] <- families[[j]]@offset(y = c(response[[j]]), w = weights)
    }
    
    ## SB: e.g. for Families 'GaussianLSS' the offsets of mu and sigma are
    ## both written into both environments of the Family for mu and sigma
    ## for the computation of ngradient()
    for (k in mods){ ## <FIXME> SB: to this double-loop outside of this j-loop?
      for (l in mods){
        if (!is.null(offset[[l]]))
          assign(names(offset)[l], families[[l]]@response(offset[[l]]),
                 environment(families[[k]]@ngradient))
      }
    }
  }
  for (j in mods){
    ## update value of nuisance parameters in families
    for (k in mods[-j]){      ## fit in comments and new call to fit must be the same 
      if (!is.null(fit[[k]])){
        ## families[[k]]@response(mboost:::fitted.mboost(fit[[k]]))
        ## <SB> set toFDboost=FALSE, to get predictions always in long format, never as matrix
        ##assign(names(fit)[k], fitted(fit[[k]], type = "response", toFDboost=FALSE),
        assign(names(fit)[k], families[[k]]@response(fitted(fit[[k]], toFDboost=FALSE)),
               environment(families[[j]]@ngradient))
      } 
    }
    ## use appropriate nu for the model
    control$nu <- nu[[j]]
    ## <FIXME> Do we need to recompute ngradient?
    ## SB: look at environment within family 
    ## ls.str(environment(families[[j]]@ngradient))
    fit[[j]] <- do.call(fun, list(formula[[names(families)[[j]]]], 
                                  timeformula = timeformula[[names(families)[[j]]]], ## <SB> timeformula for FDboost
                                  data = data, family = families[[j]],
                                  offset = "scalar", ## <SB> fixme: always use scalar offset?
                                  control=control, weights = w,
                                  ...))    
  }
  if (trace)
    do_trace(current = 1, mstart = 0,
             mstop = max(mstop),
             risk = fit[[length(fit)]]$risk())
  
  ### set up a function for iterating boosting steps
  iBoost <- function(niter) {
    start <- sapply(fit, mstop)
    mvals <- vector("list", length(niter))
    for (j in 1:length(niter)){
      mvals[[j]] <- rep(start[j] + niter[j], max(niter))
      if (niter[j] > 0)
        mvals[[j]][1:niter[j]] <- (start[j] + 1):(start[j] + niter[j])
    }
    
    ENV <- lapply(mods, function(j) environment(fit[[j]]$subset))
    for (i in 1:max(niter)){
      for (j in mods){
        ## update value of nuisance parameters
        ## use response(fitted()) as this is much quicker than 
        ## fitted(, type = response) as the latter uses predict()
        ## <SB> set toFDboost=FALSE, to get predictions always in long format, never as matrix
        for (k in mods[-j])
          ##assign(names(fit)[k], fitted(fit[[k]], type = "response", toFDboost=FALSE),
          assign(names(fit)[k], families[[k]]@response(fitted(fit[[k]], toFDboost=FALSE)),
                 environment(get("ngradient", environment(fit[[j]]$subset))))
        ## update value of u, i.e. compute ngradient with new nuisance parameters
        
        ENV[[j]][["u"]] <- ENV[[j]][["ngradient"]](ENV[[j]][["y"]], ENV[[j]][["fit"]],
                                                   ENV[[j]][["weights"]])
        # same as:
        # evalq(u <- ngradient(y, fit, weights), environment(fit[[j]]$subset))
        
        ## update j-th component to "m-th" boosting step
        fit[[j]][mvals[[j]][i]]
      }
      if (trace){
        ## which is the current risk? rev() needed to get the last
        ## list element with maximum length
        whichRisk <- names(which.max(rev(lapply(lapply(fit, function(x) x$risk()), length))))
        do_trace(current = max(sapply(mvals, function(x) x[i])),
                 mstart = ifelse(firstRun, 0, max(start)),
                 mstop = ifelse(firstRun, max(niter) + 1, max(niter)),
                 risk = fit[[whichRisk]]$risk())
      }
    }
    return(TRUE)
  }
  
  if (any(mstop > 1)){
    ## actually go for initial mstop iterations!
    firstRun <- TRUE
    tmp <- iBoost(mstop - 1)
  }
  
  firstRun <- FALSE
  
  class(fit) <- c(paste(funchar, "LSS", sep=""), "mboostLSS")
  
  ### update to a new number of boosting iterations mstop
  ### i <= mstop means less iterations than current
  ### i >  mstop needs additional computations
  ### updates take place in THIS ENVIRONMENT,
  ### some models are CHANGED!
  attr(fit, "subset") <- function(i) {
    
    i <- check(i, "mstop", names(families))
    
    msf <- mstop(fit)
    niter <- i - msf
    if (all(niter == 0)) {
      ## make nothing happen with the model
      minStart <- max(msf)
    } else {
      minStart <- min(msf[niter != 0], i[niter != 0])
    }
    
    ## check if minStart bigger than mstop of parameters that are not touched
    #if (length(msf[niter == 0]) > 0 && minStart < min(msf[niter == 0]))
    #minStart <- min(msf[niter == 0])
    
    ## reduce models first (when necessary)
    if (any(msf > minStart)){
      cf <- class(fit)
      class(fit) <- "list" ## needed to use [] operator for lists
      
      #cat("processed parameters: ", paste(names(fit[msf > minStart]),
      #                                    collapse = ", "), "\n")
      
      lapply(fit[msf > minStart],
             function(obj) obj$subset(minStart))
      
      ## remove additional boosting iterations from environments
      lapply(fit[msf > minStart], function(obj){
        evalq({xselect <- xselect[1:mstop];
        mrisk <- mrisk[1:mstop];
        ens <- ens[1:mstop];
        nuisance <- nuisance[1:mstop]},
        environment(obj$subset))
      })
      
      class(fit) <- cf
      if (trace)
        cat("Model first reduced to mstop = ", minStart, ".\n",
            "Now continue ...\n", sep ="")
    }
    
    ## now increase models (when necessary)
    if (any(i > minStart)){
      ## set negative values to 0
      ## (only applicable if some parameters do not need to be touched
      inc <- ifelse(i - minStart > 0, i - minStart, 0)
      tmp <- iBoost(inc)
    }
    
    mstop <<- i
  }
  
  ## make call in submodels nicer:
  cl <- call
  cl[[1]] <- as.name(gsub("LSS", "", cl[[1]]))
  names(cl)[names(cl) == "families"] <- "family"
  for (i in 1:length(fit)) {
    fit[[i]]$call <- cl
    ## <FIXME> This is not really nice
    fit[[i]]$call$family <- families[[i]]
  }
  
  attr(fit, "(weights)") <- weights  ## attach weights used for fitting
  
  ## update to new weights; just a fresh start
  attr(fit, "update") <- function(weights = NULL, oobweights = NULL,
                                  risk = NULL, trace = NULL, mstop = NULL) {
    if (is.null(mstop)) {
      control$mstop <- mstoparg
    } else {
      control$mstop <- mstop
    }
    if (!is.null(risk))
      control$risk <- risk
    if (!is.null(trace))
      control$trace <- trace
    ## re-use user specified offset only
    ## (since it depends on weights otherwise)
    ## this is achieved via a re-evaluation of the families argument
    FDboostLSS_fit(formula = formula, timeformula = timeformula, 
                   data = data,
                   families = eval(call[["families"]]), weights = weights,
                   control = control, fun = fun, funchar = funchar,
                   call = call, oobweights = oobweights)
  }
  attr(fit, "control") <- control
  attr(fit, "call") <- call
  attr(fit, "data") <- data
  attr(fit, "families") <- families
  return(fit)
}


#################################################################################
#' Cross-validation for FDboostLSS
#' 
#' Multidimensional cross-validated estimation of the empirical risk for hyper-parameter selection, 
#' for an object of class \code{FDboostLSS} setting the folds per default to resampling curves.  
#' 
#' @param object an object of class \code{FDboostLSS}. 
#' @param folds a weight matrix a weight matrix with number of rows equal to the number of observations. 
#' The number of columns corresponds to the number of cross-validation runs, 
#' defaults to 25 bootstrap samples, resampling whole curves  
#' @param grid a matrix of stopping parameters the empirical risk is to be evaluated for. 
#' Each row represents a parameter combination. The number of columns must be equal to the number 
#' of parameters of the GAMLSS family. Per default, make.grid(mstop(object)) is used. 
#' @param papply (parallel) apply function, defaults to \code{\link[parallel]{mclapply}}, 
#' see \code{\link[gamboostLSS]{cvrisk.mboostLSS}} for details 
#' @param trace print status information during cross-validation? Defaults to \code{TRUE}.
#' @param fun if \code{fun} is \code{NULL}, the out-of-sample risk is returned. 
#' \code{fun}, as a function of \code{object}, 
#' may extract any other characteristic of the cross-validated models. These are returned as is.
#' @param ... additional arguments passed to \code{\link[parallel]{mclapply}}.
#' 
#' @details The function \code{cvrisk.FDboostLSS} is a wrapper for 
#' \code{\link[gamboostLSS]{cvrisk.mboostLSS}} in package gamboostLSS.  
#' It overrieds the default for the folds, so that the folds are sampled on the level of curves 
#' (not on the level of single observations, which does not make sense for functional response).  
#' 
#' @return An object of class \code{cvriskLSS} (when \code{fun} was not specified), 
#' basically a matrix containing estimates of the empirical risk for a varying number 
#' of bootstrap iterations. \code{plot} and \code{print} methods are available as well as an 
#' \code{mstop} method, see \code{\link[gamboostLSS]{cvrisk.mboostLSS}}.
#' 
#' @seealso \code{\link[gamboostLSS]{cvrisk.mboostLSS}} in packge gamboostLSS. 
#' 
#' @export
## wrapper for cvrisk of gamboostLSS, specifying folds on the level of curves
cvrisk.FDboostLSS <- function(object, folds = cvLong(id = object[[1]]$id, 
                                                     weights = model.weights(object[[1]])),
                              grid = make.grid(mstop(object)),
                              papply = mclapply, trace = TRUE, 
                              fun = NULL, ...){
  
  ## message not necessary as currently only a scalar offset is possible for FDboostLSS-models
  ## if(!length(unique(object$offset)) == 1) message("The smooth offset is fixed over all folds.")
  
  class(object) <- "mboostLSS"
  
  ret <- cvrisk(object = object, folds = folds,
                grid = grid,
                papply = papply, trace = trace, 
                fun = fun, ...)
  return(ret) 
}


# #################################################################################
# #' Extract selected base-learners
# #' 
# #' Extract selected base-learners, see also \code{selected.mboostLSS}.  
# #' 
# #' @param object an object of class 'FDboostLSS' 
# #' 
# #' @return list with an entry for each distribution parameter 
# #' giving the numbers of the base-learners that were selected in each 
# #' boosting-iteration. 
# #' 
# #' @export
# selected.FDboost <- function(object, ...)
#   object$xselect()

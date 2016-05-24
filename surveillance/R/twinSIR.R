################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Function 'twinSIR' performs (penalized) maximum likelihood inference 
### for the Hoehle (2009) model. Now with REML estimation of smoothing
### parameter lambda.
###
### Copyright (C) 2008-2009 Michael Hoehle, 2008-2009, 2014 Sebastian Meyer
### $Revision: 1079 $
### $Date: 2014-10-18 01:26:00 +0200 (Sam, 18. Okt 2014) $
################################################################################

## ATTENTION: the .loglik and .score functions assume atRiskY == 1 data

######################################################################
# Log-Likelihood function
#
# PARAMS:
#  theta - parameter vector c(alpha,beta), where
#          beta also contains the baseline coefficients in the first place
#  X     - covariate matrix related to alpha, i.e. the epidemic component
#  Z     - covariate matrix related to beta, i.e. the Cox-like endemic component
#  survs - data.frame with columns id, start, stop and event
#  weights - vector of length nrow(X) indicating the number of individuals
#            with the same covariates. weights are allowed to change over time.
#            Note: it is assumed that none of the individuals covered by
#            "weights"  can have an actual event, if so they need to have their
#            own row
######################################################################

.loglik <- function(theta, X, Z, survs, weights)
{
  # Calculate epidemic (e) and endemic (h) component of the infection intensity
  eh <- .eh(theta, X, Z)
  
  # Calculate infection intensity assuming atRiskY == 1 for all rows
  lambdaNoY <- rowSums(eh)
  
  # dN Part of the loglik
  isEvent <- survs$event == 1
  events <- which(isEvent)
  intdN <- numeric(length(isEvent))   # zeros
  intdN[events] <- weights[events] * log(lambdaNoY[events])
  # here one might have got -Inf values in case of 0-intensity at an event time
  
  # lambda integral of the log-likelihood
  dt <- survs$stop - survs$start
  intlambda <- weights * lambdaNoY * dt

  # Return the log-likelihood
  loglik <- sum( intdN - intlambda )
  return(loglik)
}


######################################################################
# Penalized log-likelihood function
# Additional Params:
# lambda.smooth  - smoothing parameter
# K              - penalty matrix on the beta component
######################################################################

.ploglik <- function(theta, X, Z, survs, weights, lambda.smooth, K)
{
  loglik <- .loglik(theta, X, Z, survs, weights)
  
  if (lambda.smooth == 0) {
    return(loglik)
  }
  
  # Add penalty term and return the penalized log-likelihood
  beta <- theta[ncol(X) + seq_len(ncol(Z))]
  penalty <- lambda.smooth/2 * drop(t(beta) %*% K %*% beta)
  return(loglik - penalty)
}



######################################################################
# Score function
# Params: see .loglik
######################################################################

.score <- function(theta, X, Z, survs, weights)
{
  dimX <- dim(X)
  nRows <- dimX[1]
  px <- dimX[2]
  pz <- ncol(Z)
  isEvent <- survs$event == 1       # event indicator for the dN integral
  events <- which(isEvent)
  dt <- survs$stop - survs$start    # for the dt integral
  
  # Calculate epidemic (e) and endemic (h) component of the infection intensity
  eh <- .eh(theta, X, Z)
  h <- eh[,2,drop=TRUE]
  
  # Calculate infection intensity at event times
  lambdaEvents <- rowSums(eh[events,,drop=FALSE])
  
  score <- if (px > 0L) {
    wX <- X * weights
    part1intdN <- matrix(0, nrow = nRows, ncol = px, dimnames = dimnames(X))
    part1intdN[events,] <- wX[events,] / lambdaEvents
    part1intlambda <- wX * dt
    colSums(part1intdN - part1intlambda)
  } else NULL
  if (pz > 0L) {
    wZh <- Z * (h * weights)
    part2intdN <- matrix(0, nrow = nRows, ncol = pz, dimnames = dimnames(Z))
    part2intdN[events,] <- wZh[events,] / lambdaEvents
    part2intlambda <- wZh * dt
    part2 <- colSums(part2intdN - part2intlambda)
    score <- c(score, part2)
  }

  return(score)
}


######################################################################
# Penalized Score function
# Additional Params: see .ploglik
######################################################################

.pscore <- function(theta, X, Z, survs, weights, lambda.smooth, K, ...)
{
  score <- .score(theta, X, Z, survs, weights)
  
  if (lambda.smooth == 0) {
    return(score)
  }
  
  # Add penalty term and return the penalized Score function
  beta <- theta[ncol(X) + seq_len(ncol(Z))]
  penalty <- c(rep.int(0, ncol(X)), lambda.smooth * K %*% beta)
  return(score - penalty)
}



######################################################################
# Fisher information matrix function
# Params: see .loglik
######################################################################

.fisherinfo <- function(theta, X, Z, survs, weights)
{
  px <- ncol(X)
  pz <- ncol(Z)
  isEvent <- survs$event == 1   # event indicator
  events <- which(isEvent)
  
  # Fisher matrix calculation only incorporates data at event times!
  Xevents <- X[events,,drop = FALSE]
  Zevents <- Z[events,,drop = FALSE]
  
  # Calculate epidemic (e) and endemic (h) component of the infection intensity
  eh <- .eh(theta, Xevents, Zevents)
  h <- eh[,2,drop=TRUE]
  
  # Calculate infection intensity
  lambda <- rowSums(eh)
  
  # calculate intdN of d/dtheta log(lambda_i(t)) for all individuals with events
  wpl <- weights[events] / lambda
  dloglambda <- if (px > 0L) Xevents * wpl else NULL
  if (pz > 0L) {
    dloglambda <- cbind(dloglambda, Zevents * (h * wpl))
  }
  
  # Build the optional variation process (Martinussen & Scheike, p64)
  fisherinfo <- matrix(0, nrow=px+pz, ncol=px+pz)
  for (i in 1:nrow(dloglambda)) {
    x <- dloglambda[i,,drop=FALSE]  # single-ROW matrix
    fisherinfo <- fisherinfo + crossprod(x) # t(x) %*% x
  }
  
  return(fisherinfo) 
}


######################################################################
# Fisher information matrix function
# Additional Params: see .ploglik
######################################################################

.pfisherinfo <- function(theta, X, Z, survs, weights, lambda.smooth, K)
{
  fisherinfo <- .fisherinfo(theta, X, Z, survs, weights)
  
  if (lambda.smooth == 0) {
    return(fisherinfo)
  }
  
  # Add penalty term and return the penalized Fisher information matrix
  penalty <- matrix(0, ncol=ncol(fisherinfo), nrow=nrow(fisherinfo))
  zIndex <- ncol(X) + seq_len(ncol(Z))
  penalty[zIndex,zIndex] <- lambda.smooth * K
  return(fisherinfo + penalty)
}

######################################################################
# Marginal likelihood of the log(smoothing) parameter as given
# by a Laplace approximation c.f. Kneib & Fahrmeir (2006), p.9.
# or Cai et al (2002)
#
# Params:
#  log.lambda.smooth - log parametrization to ensure positive value of
#                      lambda.smooth
#  theta - fixed regression parameters
#  X - design matrix of additive part
#  Z - design matrix of multiplicative part
#  survs - the data.frame containing the data in survs format
#  weights - for weighting individual entries
#  K       - smoother matrix
#
# Returns:
#  value of lmarg
######################################################################
 
.lmarg.lambda <- function(log.lambda.smooth, theta, X, Z, survs, weights, K) {
  #Contribution of the penalized likelihood
  loglik <- .ploglik(theta, X, Z, survs, weights, exp(log.lambda.smooth), K)
  #Laplace approximation using TP representation
  H <- .pfisherinfo(theta, X, Z, survs, weights, exp(log.lambda.smooth), K)
  beta <- theta[ncol(X) + seq_len(ncol(Z))]
  #[Q]: Extract baseline terms from model and translate into
  #TP-spline setting, i.e. a B-spline of 0th order is assumed
  baselineIdx <- grep("cox\\(logbaseline.*\\)",dimnames(Z)[[2]])
  b <- diff(beta[baselineIdx])
  laplace <- 1/2*(length(b)-1)*log.lambda.smooth - 1/2*log(det(H))
  return(loglik + laplace)
}


######################################################################
# Model fitter. Prepares everything and uses optim's (L-)BFGS(-B) to
# maximize the (penalized) log-likelihood.
######################################################################

twinSIR <- function (formula, data, weights, subset,
    knots = NULL, nIntervals = 1, lambda.smooth = 0, penalty = 1,
    optim.args = list(), model = TRUE, keep.data = FALSE)
{
  cl <- match.call()
  
  ## Verify that 'data' inherits from "epidata"
  data <- eval(cl$data, parent.frame())
  if (!inherits(data, "epidata")) {
    stop("'data' must inherit from class \"epidata\"")
  }
  
  ## Extract the time range of the epidemic
  timeRange <- attr(data, "timeRange")
  minTime <- timeRange[1L]
  maxTime <- timeRange[2L]
  
#   ## NOTE: modification of 'data' has no effect with the current evaluation
#   ##       of model.frame in the parent.frame() as the original 'data' will 
#   ##       be used.
#   ## Impute blocks for 'knots', which are not existing stop times
#   if (is.vector(knots, mode = "numeric")) {
#     insideKnot <- (knots > minTime) & (knots < maxTime)
#     if (any(!insideKnot)) {
#       warning("only 'knots' inside the observation period are considered")
#     }
#     knots <- sort(knots[insideKnot])
#     data <- intersperse(data, knots)
#   }
  
  
  ############################
  ### Build up model.frame ### (this is derived from the coxph function)
  ############################
  
  mfnames <- c("", "formula", "data", "weights", "subset")
  mf <- cl[match(mfnames, names(cl), nomatch = 0L)]
  mf$id <- as.name("id")
  mf$atRiskY <- as.name("atRiskY")
  mf$subset <- if (is.null(mf$subset)) {
      call("==", mf$atRiskY, 1)
    } else {
      call("&", mf$subset, call("==", mf$atRiskY, 1))
    }
  if(length(formula) == 2L) { # i.e. no response specified
    formula[3L] <- formula[2L]
    formula[[2L]] <- quote(cbind(start, stop, event))
  }
  mf$na.action <- as.name("na.fail")
  special <- c("cox")
  Terms <- terms(formula, specials = special, data = data, keep.order = FALSE)
  mf$formula <- Terms
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  
  ###########################################################
  ### Check arguments and extract components of the model ###
  ###########################################################
  
  ## Extract and check 'weights'
  weights <- model.extract(mf, "weights")
  if (is.null(weights)) {
    weights <- rep(1, nrow(mf))
    names(weights) <- attr(mf, "row.names")
  } else {
    if (!is.vector(weights, mode="numeric")) {
      stop("'weights' must be a numeric vector")
    }
    if (any(weights < 0)) {
      stop("negative 'weights' not allowed")
    }
  }
  
  ## Extract the response
  response <- model.response(mf)
  survs <- data.frame(id = model.extract(mf, "id"), start = response[,1L],
                      stop = response[,2L], event = response[,3L],
                      check.names = FALSE, stringsAsFactors = FALSE)
  attr(survs, "eventTimes") <- attr(data, "eventTimes")
  attr(survs, "timeRange") <- timeRange
  
  ## Check specified baseline intervals
  if (is.null(knots) && isScalar(nIntervals)) {
    knots <- if (nIntervals == 1) {
      numeric(0)
    } else if (nIntervals > 1) {
      quantile(attr(survs, "eventTimes"),
        probs = seq(from=0, to=1, length.out=nIntervals+1)[-c(1,nIntervals+1)],
        type = 1, names = FALSE)
    } else {
      stop("'nIntervals' must be a single number >= 1")
    }
  } else if (is.vector(knots, mode = "numeric")) {
    isInsideKnot <- (knots > minTime) & (knots < maxTime)
    if (any(!isInsideKnot)) {
      warning("only 'knots' inside the observation period are considered")
      knots <- knots[isInsideKnot]
    }
    isStopKnot <- knots %in% unique(survs$stop)
    if (any(!isStopKnot)) {
      stop("ATM, 'knots' must be a subset of the 'stop' times where at least ",
           "one individual is susceptible")
#       nNonStopKnots <- sum(!isStopKnot)
#       warning(
#         sprintf(ngettext(nNonStopKnots,
#           paste("%d of 'knots' has been ignored due to no susceptible",
#                 "individuals at this time point: "),
#           paste("%d 'knots' have been ignored due to no susceptible",
#                 "individuals at those time points: ")), nNonStopKnots),
#         knots[!stopKnot]
#       )
#       knots <- knots[stopKnot]
    }
    knots <- sort(knots)
  } else {
    stop("'knots' (a numeric vector) or 'nIntervals' (a single number) ",
         "must be specified")
  }
  
  intervals <- c(minTime, knots, maxTime)
  nIntervals <- length(intervals) - 1L
  message(
    sprintf(ngettext(nIntervals,
                     "Initialized %d log-baseline interval:  ",
                     "Initialized %d log-baseline intervals:  "),
            nIntervals), 
    paste(format(intervals, trim = TRUE), collapse=" ")
  )
  
  ## Extract the two parts of the design matrix:
  ## Z contains the Cox part, X contains the epidemic part, there's no intercept
  des <- read.design(mf, Terms)
  X <- des$X; px <- ncol(X)
  Z <- des$Z
  
  ## Add variables for the piecewise constant baseline to Z (if requested)
  if (nIntervals == 1L) {
    nEvents <- length(attr(survs, "eventTimes"))
    if (attr(Terms, "intercept") == 1) Z <- cbind("cox(logbaseline)" = 1, Z)
  } else {   # we have more than one baseline interval/parameter
    intervalIndices <- findInterval(survs$start, intervals,
                                    rightmost.closed = FALSE)
    intervalNumbers <- seq_len(nIntervals)
    baselineVars <- sapply(intervalNumbers, function(i) intervalIndices == i)
    dimnames(baselineVars) <- list(NULL,
        paste("cox(logbaseline.", intervalNumbers, ")", sep=""))
    Z <- cbind(baselineVars, Z)
    nEvents <- as.vector(table(factor(intervalIndices[survs$event == 1],
                                      levels = seq_len(nIntervals))))
  }
  pz <- ncol(Z)
  
  ## Check that we have at least one parameter
  if (pz == 0L && px == 0L) {
    stop("nothing to do: neither a baseline nor covariates have been specified")
  }
  
  ## Check lambda.smooth
  if (!isScalar(lambda.smooth)) {
    stop("'lambda.smooth' must be scalar")
  }
  if (lambda.smooth != 0 && pz == 0L) {
    lambda.smooth <- 0
    message("Note: 'lambda.smooth' was set to 0, because there was no endemic ",
            "component in the formula.")
  }

  ## Setup penalty matrix
  if (isScalar(penalty)) {
    K <- matrix(0, ncol = pz, nrow = pz)
    if (lambda.smooth != 0 && nIntervals > 1L) {
      # do we have equidistant knots?
      knotSpacings <- diff(intervals)
      #equidistant <- all(sapply(knotSpacings[-1], function(x) isTRUE(all.equal(x,knotSpacings[1]))))
      equidistant <- isTRUE(all.equal(diff(knotSpacings), rep.int(0,nIntervals-1)))
      if (equidistant) { # K = D'D only works for equidistant knots
        # difference matrix of order 'penalty'
        D <- diff(diag(nIntervals), differences=penalty)
        K[intervalNumbers,intervalNumbers] <- crossprod(D)   # t(D) %*% D
      } else { # special weighting scheme for the non-equidistant case
        if (penalty != 1) {
          stop("ATM, non-equidistant knots only work for 1st order penalty")
        }
        #Use Fahrmeir & Lang (2001), p.206
        invdelta <- 1/diff(intervals) * mean(diff(intervals))
        #Use Fahrmeir & Lang (2001), p.206
        for (i in 1:(nIntervals)) {
          idx2 <- cbind(j=c(-1,1) + i, deltaidx=i+c(-1,0),fac=c(-1,-1))
          idx2 <- idx2[idx2[,"j"] > 0 & idx2[,"j"] <= nIntervals,,drop=FALSE]
          #Off diagonal elements
          K[i, idx2[,"j"]] <- invdelta[idx2[,"deltaidx"]] * idx2[,"fac"]
          #Diagonal element
          K[i, i] <- sum(invdelta[idx2[,"deltaidx"]])
        }
        message("Note: non-equidistant knots. Using penalization matrix ",
                "correcting for distance between knots.\n")
#        print(K)
#        browser()
      }
    }
  } else if (is.matrix(penalty) && ncol(penalty) == pz && nrow(penalty) == pz) {
    K <- penalty
  } else {
    stop("'penalty' must either be a single number or a square matrix of ",
         "dimension ", pz, "x", pz, ", fitting the number of unknown ",
         "parameters in the endemic component (baseline and covariates)")
  }
  
  ## Check that optim.args is a list
  if (!is.list(optim.args)) {
    stop("'optim.args' must be a list")
  }
  
  ## Check start value for theta
  if (!is.null(optim.args[["par"]])) {
    if (!is.vector(optim.args$par, mode="numeric")) {
      stop("'optim.args$par' must be a numeric vector or NULL")
    }
    if (length(optim.args$par) != px + pz) {
      stop(gettextf(paste("'optim.args$par' (%d) does not have the same length",
                          "as the number of unknown parameters (%d + %d = %d)"),
                    length(optim.args$par), px, pz, px + pz))
    }
  } else {
    optim.args$par <- c(rep.int(1, px), rep.int(0, pz))
  }
  message("Initial parameter vector:  ", paste(optim.args$par, collapse=" "))
  
  ## Set names for theta
  names(optim.args$par) <- c(colnames(X), colnames(Z))


  ####################
  ### Optimization ###
  ####################
  
  ## Configuring the optim procedure (check optim.args)
  optimControl <- list(trace = 1, fnscale = -1, maxit = 300, factr = 1e7)
  optimControl[names(optim.args[["control"]])] <- optim.args[["control"]]
  optim.args$control <- optimControl
  optimArgs <- list(par = optim.args$par, fn = .ploglik, gr = .pscore,
                    X = X, Z = Z, survs = survs, weights = weights,
                    lambda.smooth = lambda.smooth, K = K, method = "L-BFGS-B",
                    lower = c(rep(0,px), rep(-Inf,pz)), upper = rep(Inf,px+pz),
                    control = list(), hessian = FALSE)
  namesOptimArgs <- names(optimArgs)
  namesOptimUser <- names(optim.args)
  optimValid <- namesOptimUser %in% namesOptimArgs
  optimArgs[namesOptimUser[optimValid]] <- optim.args[optimValid]
  if (any(!optimValid))
    warning("unknown names in optim.args: ",
      paste(namesOptimUser[!optimValid], collapse = ", "))
  if (! "method" %in% namesOptimUser && px == 0L) {
    optimArgs$method <- "BFGS"
  }
  if (optimArgs$method != "L-BFGS-B") {
    optimArgs$lower <- -Inf
    optimArgs$upper <- Inf
  }

  #Fit model using fixed smoothing parameter or use mixed model
  #representation to estimate lambda.smooth using marginal likelihood
  if (lambda.smooth == -1) {
    if (isScalar(penalty) && penalty == 1) {
      ###################################################################
      ##TODO: Need to check for B-spline (?). Move options into ctrl obj
      ###################################################################
      #Iterative procedure where we change between optimizing regression
      #parameters given fixed smoothing parameter and optimizing the
      #smoothing parameter given fixed regression parameters (Gauss-Seidel)
      #procedure. The tuning parameters (5) could go into the control object.
      lambda.smooth <- 5
      reltol <- 1e-2
      maxit <- 25

      #Parameters for keeping track of the iterations
      lambda.smoothOld <- 1e99
      iter <- 0

      #Loop until relative convergence or max-iteration reached
      while ((abs(lambda.smooth-lambda.smoothOld)/lambda.smoothOld > reltol) &
             (iter < maxit)) {
        #Iteration begins
        iter <- iter + 1
        if (optimControl$trace > 0) {
          cat("==> Iteration ",iter," of Gauss-Seidel maximization. lambda.smooth = ",lambda.smooth,"\n")
        }
        
        #Step 1 - maximize (alpha,beta) with fixed lambda
        optimArgs$lambda.smooth <- lambda.smooth
        optimRes <- do.call("optim", optimArgs)
        theta <- optimRes$par
        optimArgs$par <- theta #better start value the next time
        
        #Step 2 - maximize log(lambda) with fixed (alpha,beta)
        optimLambda <- optim(log(lambda.smooth), .lmarg.lambda, control=list(fnscale=-1,trace=1),method="BFGS",
                             theta=theta, X=X, Z=Z, survs=survs, weights=weights, K=K)
        lambda.smoothOld <- lambda.smooth
        lambda.smooth <- exp(optimLambda$par)
      }
      #Done, update optimArgs with new smoothing parameter
      optimArgs$lambda.smooth <- lambda.smooth
    } else {
      stop("REML estimation using TP-splines only works for 1st order differences.")
    }
  }

  ## Call optim with the arguments above (including the news smoothing param)
  optimRes <- do.call("optim", optimArgs)

  ##############
  ### Return ###
  ##############

  ## Set up list object to be returned
  fit <- list(coefficients = optimRes$par, lambda.smooth = lambda.smooth, loglik = optimRes$value,
              counts = optimRes$counts, converged = (optimRes$convergence == 0))
              
  ## If requested, add observed fisher info (= negative hessian at maximum)
  if (!is.null(optimRes$hessian)) {
    fit$fisherinfo.observed <- -optimRes$hessian
  }
  
  ## Add own (exact) fisher info computation
  fit$fisherinfo <- .pfisherinfo(theta = fit$coefficients, X = X, Z = Z,
                                 survs = survs, weights = weights,
                                 lambda.smooth = lambda.smooth, K = K)
  
  ## Add 'method'
  fit$method <- optimArgs$method
  
  ## Append further information
  fit$intervals <- intervals
  fit$nEvents <- nEvents
  if (model) {
    fit$model <- list(
      survs = survs, X = X, Z = Z, weights = weights,
      lambda.smooth = lambda.smooth, K = K,
      f = attr(data, "f")[match(colnames(X), names(attr(data, "f")), nomatch=0)],
      w = attr(data, "w")[match(colnames(X), names(attr(data, "w")), nomatch=0)]
    )
  }
  if (keep.data) {
    fit$data <- data
  }
  fit$call <- cl
  fit$formula <- formula(Terms)
  fit$terms <- Terms
  
  ## Return object of class "twinSIR"
  class(fit) <- "twinSIR"
  return(fit)
}

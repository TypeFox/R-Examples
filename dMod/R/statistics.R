#' Profile-likelihood (PL) computation
#' 
#' @param obj Objective function \code{obj(pars, fixed, ...)} returning a list with "value",
#' "gradient" and "hessian".
#' @param pars Parameter vector corresponding to the log-liklihood optimum.
#' @param whichPar Numeric or character. The parameter for which the profile is computed.
#' @param alpha Numeric, the significance level based on the chisquare distribution with df=1
#' @param limits Numeric vector of length 2, the lower and upper deviance from the original 
#' value of \code{pars[whichPar]}
#' @param stepControl List of arguments controlling the step adaption. Defaults to 
#' \code{list(stepsize = 1e-4, min = 0, max = Inf, atol = 1e-1, rtol = 1e-1, limit = 100)}
#' @param algoControl List of arguments controlling the fast PL algorithm. defaults to
#' \code{list(gamma = 1, W = c("hessian", "identity"), reoptimize = FALSE, correction = 1, reg = 1e-6)}
#' @param optControl List of arguments controlling the \code{trust()} optimizer. Defaults to
#' \code{list(rinit = .1, rmax = 10, iterlim = 10, fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))}.
#' See \link{trust} for more details.
#' @param verbose Logical, print verbose messages.
#' @param ... Arguments going to obj()
#' @details Computation of the profile likelihood is based on the method of Lagrangian multipliers
#' and Euler integration of the corresponding differential equation of the profile likelihood paths.
#' 
#' \code{algoControl}: Since the Hessian which is needed for the differential equation is frequently misspecified, 
#' the error in integration needs to be compensated by a correction factor \code{gamma}. Instead of the
#' Hessian, an identity matrix can be used. To guarantee that the profile likelihood path stays on
#' the true path, each point proposed by the differential equation can be used as starting point for
#' an optimization run when \code{reoptimize = TRUE}. The correction factor \code{gamma} is adapted
#' based on the amount of actual correction. If this exceeds the value \code{correction}, \code{gamma} is
#' reduced. In some cases, the Hessian becomes singular. This leads to problems when inverting the
#' Hessian. To avoid this problem, the matrix \code{reg*Id} is added to the Hessian.
#' 
#' \code{stepControl}: The Euler integration starts with \code{stepsize}. In each step the predicted change
#' of the objective function is compared with the actual change. If this is larger than \code{atol}, the
#' stepsize is reduced. For small deviations, either compared the abolute tolerance \code{atol} or the
#' relative tolerance \code{rtol}, the stepsize may be increased. \code{max} and \code{min} are upper and lower
#' bounds for \code{stepsize}. \code{limit} is the maximum number of steps that are take for the profile computation.
#' 
#' @return Named list of length one. The name is the parameter name. The list enty is a
#' matrix with columns "value" (the objective value), "constraint" (deviation of the profiled paramter from
#' the original value), "stepsize" (the stepsize take for the iteration), "gamma" (the gamma value employed for the
#' iteration), one column per parameter (the profile paths).
#' @examples 
#' 
#' \dontrun{
#' ## ----------------------
#' ## Example 1 
#' ## ----------------------
#' trafo <- c(a = "exp(loga)", b = "exp(logb)",c = "exp(loga)*exp(logb)*exp(logc)")
#' p <- P(trafo) 
#' obj <- function(pOuter, fixed = NULL) 
#'    constraintL2(p(pOuter, fixed), c(a =.1, b = 1, c = 10), 1)
#'     
#' ini <- c(loga = 1, logb = 1, logc = 1)   
#' myfit <- trust(obj, ini, rinit=1, rmax=10)   
#' profiles <- sapply(1:3, function(i) 
#'    profile(obj, myfit$argument, whichPar = i, limits = c(-5, 5), 
#'                  algoControl=list(gamma=1, reoptimize=FALSE), verbose=TRUE))
#' plotProfile(profiles)
#' plotPaths(profiles)
#' 
#' ## ----------------------------
#' ## Example 2
#' ## ----------------------------
#' trafo <- c(a = "exp(loga)", b = "exp(logb)",c = "exp(loga)*exp(logb)*exp(logc)")
#' p <- P(trafo)
#' obj <- function(pOuter, fixed = NULL, sigma) 
#'   constraintL2(p(pOuter, fixed), c(a =.1, b = 1, c = 10), 1) +
#'   constraintL2(pOuter, mu = c(loga = 0, logb = 0), sigma = sigma, fixed = fixed)
#' 
#' 
#' ini <- c(loga = 1, logb = 1, logc = 1)
#' myfit <- trust(obj, ini[-1], rinit=1, rmax=10, fixed = ini[1], sigma = 10)
#' profiles.approx <- sapply(1:2, function(i) 
#'   profile(obj, myfit$argument, whichPar = i, limits = c(-10, 10), 
#'                 algoControl=list(gamma=1, reoptimize=FALSE), 
#'                 verbose=TRUE, fixed = ini[1], sigma = 10))
#' profiles.exact  <- sapply(1:2, function(i) 
#'   profile(obj, myfit$argument, whichPar = i, limits = c(-10, 10), 
#'                 algoControl=list(gamma=0, reoptimize=TRUE), 
#'                 verbose=TRUE, fixed = ini[1], sigma = 10))
#' 
#' plotProfile(profiles.approx, profiles.exact)
#' }
#' @export
#' @import trust
profile <- function(obj, pars, whichPar, alpha = 0.05, 
                          limits = c(lower = -Inf, upper = Inf), 
                          stepControl = NULL, 
                          algoControl = NULL,
                          optControl  = NULL,
                          verbose = FALSE,
                          ...) {
  
  sControl <- list(stepsize = 1e-4, min = 0, max = Inf, atol = 1e-1, rtol = 1e-1, limit = 100)
  aControl <- list(gamma = 1, W = c("hessian", "identity"), reoptimize = FALSE, correction = 1, reg = 1e-6)
  oControl <- list(rinit = .1, rmax = 10, iterlim = 10, fterm = sqrt(.Machine$double.eps), mterm = sqrt(.Machine$double.eps))
  
  if(!is.null(stepControl)) sControl[match(names(stepControl), names(sControl))] <- stepControl
  if(!is.null(algoControl)) aControl[match(names(algoControl), names(aControl))] <- algoControl
  if(!is.null(optControl )) oControl[match(names(optControl), names(oControl ))] <- optControl
    
  
  if(is.character(whichPar)) whichPar <- which(names(pars) == whichPar)
  whichPar.name <- names(pars)[whichPar]
  if(any(names(list(...)) == "fixed")) fixed <- list(...)$fixed else fixed <- NULL
  
  
  ## Functions needed during profile computation -----------------------
  obj.opt <- obj
  obj.prof <- function(p, ...) {
    out <- obj(p, ...)
    # Substitute hessian by the identity matrix
    Id <- diag(1, length(out$gradient))
    colnames(Id) <- rownames(Id) <- names(out$gradient)
    W <- match.arg(aControl$W[1], c("hessian", "identity"))
    out$hessian <- switch(W,
                          "hessian" = out$hessian,
                          "identity" = Id)
    return(out)    
  }
  constraint <- function(p) {
    value <- p[whichPar] - pars[whichPar]
    gradient <- rep(0, length(p))
    gradient[whichPar] <- 1
    return(list(value = value, gradient = gradient))
  }
  lagrange <- function(y) {
    
    # initialize values
    p <- y
    lambda <- 0
    out <- obj.prof(p, ...)
    g.original <- constraint(p)
    
    # evaluate derivatives and constraints
    g     <- direction * g.original$value
    gdot  <- direction * g.original$gradient
    ldot  <- out$gradient
    lddot <- out$hessian  + diag(aControl$reg, length(out$gradient)) 
    
    # compute rhs of profile ODE
    M <- rbind(cbind(lddot, gdot), 
               matrix(c(gdot, 0), nrow=1))
    
    v <- c(-rep(gamma, length(p))*(ldot + lambda*gdot), 1)
    v0 <- c(-rep(0, length(p))*(ldot + lambda*gdot), 1)
    
    dy <- try(as.vector(solve(M)%*%v)[1:length(p)], silent=FALSE)
    dy0 <- try(as.vector(solve(M)%*%v0)[1:length(p)], silent=FALSE)
    
    if(!inherits(dy, "try-error")) {
      names(dy) <- names(y) 
      correction <- sqrt(sum((dy-dy0)^2))/sqrt(sum(dy^2))
    } else {
      dy <- NA
      correction <- 0
      warning("Impossible to invert Hessian. Trying to optimize instead.")
    }
    
    return(list(dy = dy, value = out$value, gradient = out$gradient, correction = correction))
    
  }
  doIteration <- function() {
    
    optimize <- aControl$reoptimize
    # Check for error in evaluation of lagrange()
    if(is.na(dy[1])) {
      #cat("Evaluation of lagrange() not successful. Will optimize instead.\n")
      optimize <- TRUE
      y.try <- y
      y.try[whichPar] <- y[whichPar] + direction*stepsize
      rinit <- oControl$rinit
    } else {
      dy.norm <- sqrt(sum(dy^2))
      rinit <- min(c(oControl$rinit, 3*dy.norm))
      y.try <- y + dy
    }
    
    # Do reoptimization if requested or necessary
    if(optimize) {      
      parinit.opt <- y.try[-whichPar]
      fixed.opt <- c(fixed, y.try[whichPar])
          
      arglist <- c(list(objfun = obj.opt, parinit = parinit.opt, fixed = fixed.opt, rinit = rinit), 
                   oControl[names(oControl)!="rinit"],
                   list(...)[names(list(...)) != "fixed"])
      
      myfit <- try(do.call(trust::trust, arglist), silent=FALSE)
      if(!inherits(myfit, "try-error")) {
        y.try[names(myfit$argument)] <- as.vector(myfit$argument)  
      } else {
        warning("Optimization not successful. Profile may be erroneous.")
      }
      
    }
    
    return(y.try)
    
  }
  doAdaption <- function() {
    
    lagrange.out.try <- lagrange(y.try)
    
    # Predicted change of the objective value
    dobj.pred <- sum(lagrange.out$gradient*(y.try-y))
    dobj.fact <- lagrange.out.try$value - lagrange.out$value
    correction <- lagrange.out.try$correction
    
    # Gamma adaption based on amount of actual correction
    if(correction > aControl$correction) gamma <- gamma/2
    if(correction < 0.5*aControl$correction) gamma <- min(c(aControl$gamma, gamma*2))
    
    # Stepsize adaption based on difference in predicted change of objective value
    if(abs(dobj.fact - dobj.pred) > sControl$atol) {
      stepsize <- max(c(stepsize/1.5, sControl$min))
    }
    if(abs(dobj.fact - dobj.pred) < .3*sControl$atol | abs((dobj.fact - dobj.pred)/dobj.fact) < .3*sControl$rtol) {
      stepsize <- min(c(stepsize*2, sControl$max))
    }
    
    # Compute progres
    diff.thres <- diff.steps <- diff.limit <- 0
    if(threshold < Inf)
      diff.thres <- 1 - max(c(0, min(c(1, (threshold - lagrange.out.try$value)/delta))))
    if(sControl$limit < Inf)
      diff.steps <- i/sControl$limit
    diff.limit <- switch(as.character(sign(constraint.out$value)),
                         "1"  = 1 - (limits[2] - constraint.out$value)/limits[2],
                         "-1" = diff.limit <- 1 - (limits[1] - constraint.out$value)/limits[1],
                         "0"  = 0)
    
    percentage <- max(c(diff.thres, diff.steps, diff.limit))*100
    progressBar(percentage)
    
    ## Verbose
    if(verbose) {
      #cat("diff.thres:", diff.thres, "diff.steps:", diff.steps, "diff.limit:", diff.limit)
      myvalue <- format(substr(lagrange.out$value  , 0, 8), width=8)
      myconst <- format(substr(constraint.out$value, 0, 8), width=8)
      mygamma <- format(substr(gamma               , 0, 8), width=8)
      cat("\tvalue:", myvalue, "constraint:", myconst, "gamma:", mygamma) 
    }
    
    
    
    return(list(lagrange = lagrange.out.try, stepsize = stepsize, gamma = gamma))
    
    
  }
  
  ## Compute profile -------------------------------------------------
    
  # Initialize profile
  direction <- 1
  gamma <- aControl$gamma
  stepsize <- sControl$stepsize
  ini <- pars
  
  lagrange.out <- lagrange(ini)
  constraint.out <- constraint(pars)
  
  delta <- qchisq(1-alpha, 1)
  threshold <- lagrange.out$value + delta
  
  out <- c(value = lagrange.out$value, constraint = as.vector(constraint.out$value), stepsize = stepsize, gamma = gamma, ini)
    
  # Compute right profile
  cat("Computer right profile\n")
  direction <- 1
  gamma <- aControl$gamma
  stepsize <- sControl$stepsize
  y <- ini
  
  lagrange.out <- lagrange.out
  constraint.out <- constraint.out
 
  i <- 0 
  while(i < sControl$limit) {
    
    ## Iteration step
    dy <- stepsize*lagrange.out$dy
    y.try <- doIteration()
    out.try <- doAdaption()
    
    ## Set values
    y <- y.try
    lagrange.out <- out.try$lagrange
    constraint.out <- constraint(y.try)
    stepsize <- out.try$stepsize
    gamma <- out.try$gamma
    
    ## Return values 
    out <- rbind(out, 
                 c(value = lagrange.out$value, constraint = as.vector(constraint.out$value), stepsize = stepsize, gamma = gamma, y))
    
    
    if(lagrange.out$value > threshold | constraint.out$value > limits[2]) break
    
    i <- i + 1
    
    
  }
  
  # Compute left profile
  cat("\nComputer left profile\n")
  direction <- -1
  gamma <- aControl$gamma
  stepsize <- sControl$stepsize
  y <- ini
  
  lagrange.out <- lagrange(ini)
  constraint.out <- constraint(pars)
    
  i <- 0
  while(i < sControl$limit) {
    
    ## Iteration step
    dy <- stepsize*lagrange.out$dy
    y.try <- doIteration()
    out.try <- doAdaption()
    
    ## Set values
    y <- y.try
    lagrange.out <- out.try$lagrange
    constraint.out <- constraint(y.try)
    stepsize <- out.try$stepsize
    gamma <- out.try$gamma
    
    
    ## Return values
    out <- rbind(c(value = lagrange.out$value, constraint = as.vector(constraint.out$value), stepsize = stepsize, gamma = gamma, y), 
                 out)
    
    if(lagrange.out$value > threshold  | constraint.out$value < limits[1]) break
    
    i <- i + 1
    
  }
  
  out.list <- list(out)
  names(out.list) <- whichPar.name
  
  return(out.list)
  

}

#' Progress bar
#' 
#' @param percentage Numeric between 0 and 100
#' @param size Integer, the size of the bar print-out
#' @param number Logical, Indicates whether the percentage should be printed out.
#' @export
progressBar <- function(percentage, size = 50, number = TRUE) {
  
  if(percentage < 0) percentage <- 0
  if(percentage > 100) percentage <- 100
  
  out <- paste("\r|", paste(rep("=", round(size*percentage/100)), collapse=""), paste(rep(" ", size-round(size*percentage/100)), collapse=""), "|", sep="")
  cat(out)
  if(number) cat(format(paste(" ", round(percentage), "%", sep=""), width=5))
  
}


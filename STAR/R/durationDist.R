invgaussMLE <- function (yi, ni = numeric(length(yi)) + 1,
                         si = numeric(length(yi)) + 1,
                         parameterization = "sigma2"
                         ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  if (any(yi < 0)) 
    stop("yi elements must be non-negative")
  yi <- as.numeric(yi)
  if (!identical(length(yi), length(ni))) 
    stop("yi and ni should have the same length")
  if (any(ni < 0)) 
    stop("ni elements must be non-negative")
  if (!identical(class(ni), "integer") && !identical(ni, round(ni))) 
    stop("ni should be a vector of positive integers")
  if (!identical(length(yi), length(si))) 
    stop("yi and si should have the same length")
  if (any(si < 0)) 
    stop("si elements must be non-negative")
  if (!identical(si, round(si))) 
    stop("si should be a vector of positive integers")
  if (any(si > ni)) 
    stop("si elements should not be greater than ni elements")
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n", 
                   "  component 1 is the log of the mu parameter,\n", 
                   "  component 2 is the log of the sigma2 parameter,\n", 
                   "using the 'sigma2' parameterization of the inverse Gaussian.\n")
      cat(txt)
    } else {
      mu <- exp(p[1])
      sigma2 <- exp(p[2])
      -(ifelse(s.dot > 0,
               sum(dinvgauss(yi[si > 0], mu, sigma2, log = TRUE) * si[si > 0], 0)
               ) +
        ifelse(c.dot > 0,
               sum(pinvgauss(yi[ci > 0], mu, sigma2, lower.tail = FALSE,log.p = TRUE) * ci[ci > 0]), 0))
    }
  }
  s.dot <- sum(si)
  n.dot <- sum(ni)
  ci <- ni - si
  c.dot <- sum(ci)
  if (s.dot == n.dot) {
    mu.hat <- weighted.mean(yi, ni)
    inv.mean <- weighted.mean(1/yi, ni)
    sigma2.hat <- inv.mean - 1/mu.hat
    estimate <- c(mu.hat, sigma2.hat)
    observedI <- matrix(c(n.dot/(sigma2.hat * mu.hat^3), 
                          0, 0, n.dot/sigma2.hat^2/2), nrow = 2, byrow = TRUE)
    se <- sqrt(diag(solve(observedI)))
    l <- -minusLogLik(log(estimate))
  } else {
    if (s.dot >= 10) {
      mu.hat <- weighted.mean(yi, si)
      inv.mean <- weighted.mean(1/yi, si)
      sigma2.hat <- inv.mean - 1/mu.hat
    } else {
      mu.hat <- weighted.mean(yi, ni)
      inv.mean <- weighted.mean(1/yi, ni)
      sigma2.hat <- inv.mean - 1/mu.hat
    }
    mleFit <- optim(par=log(c(mu.hat, sigma2.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian = TRUE)
    estimate <- exp(mleFit$par)
    newVar <- (1/estimate) %o% (1/estimate)
    observedI <- mleFit$hessian * newVar
    se <- sqrt(diag(solve(observedI)))
    l <- -mleFit$value
  }
  if (parameterization == "sigma2") {
    names(estimate) <- c("mu", "sigma2")
    names(se) <- c("mu", "sigma2")
    rFct <- function(mu, sigma2) -minusLogLik(log(c(mu, sigma2))) - l
  } else {
    boundary.hat <- (sigma2.hat)^(-0.5)
    mu.hat <- mu.hat/boundary.hat
    if (s.dot == n.dot) {
      estimate <- c(mu.hat, boundary.hat)
      observedI <- observedI * matrix(c(boundary.hat^2, 
                                        -2/boundary.hat^2, -2/boundary.hat^2, 4/boundary.hat^6), 
                                      nrow = 2, byrow = TRUE)
      se <- se * c(1/boundary.hat, boundary.hat^3/2)
    } else {
      estimate <- c(estimate[1] * sqrt(estimate[2]), 1/sqrt(estimate[2]))
      newVar <- newVar * (c(estimate[2], -2/estimate[2]^3) %o% 
                          c(estimate[2], -2/estimate[2]^3))
      observedI <- mleFit$hessian * newVar
      se <- sqrt(diag(solve(observedI)))
    }
    names(estimate) <- c("mu", "boundary")
    names(se) <- c("mu", "boundary")
    rFct <- function(mu, boundary) -minusLogLik(log(c(mu * boundary, 1/boundary^2))) - l
  }
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct, 
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}

coef.durationFit <- function(object,...) object$estimate
logLik.durationFit <- function(object,...) object$logLik

gammaMLE <- function(yi,
                     ni = numeric(length(yi))+1,
                     si = numeric(length(yi))+1,
                     scale = TRUE
                     ) {
  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi < 0)) stop("yi elements must be non-negative")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")
  
  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the shape parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
                   )
      cat(txt)
    } else {
      shape <- exp(p[1])
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dgamma(yi[si>0],shape=shape,scale=scale,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(pgamma(yi[ci>0],shape=shape,
                                  scale=scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  ## Define function logProfileLik returning the log of the
  ## profiled likelihood for the shape parameter
  logProfileLik <- function(shape) {
    if (shape <= 0) return(-Inf)
    term1 <- s.dot*shape*(log(shape) - log(scale.hat) -1)
    term2 <- (shape - 1) * sum(ni * log(yi))
    term3 <- - s.dot * lgamma(shape)
    return(term1 + term2 + term3)
  }
  
  if (s.dot == n.dot) {
    ## no censored event
    ## Get the empirical mean which is the MLE of parameter scale
    scale.hat <- weighted.mean(yi, ni)
    s2 <- weighted.mean(yi^2, ni) - scale.hat^2
    shape.hat <- scale.hat^2 / s2
    shape.lower <- floor(shape.hat)
    shape.upper <- ceiling(shape.hat)
    shape.hat <- optimize(logProfileLik,
                          lower = shape.lower,
                          upper = shape.upper,
                          maximum = TRUE)$maximum
    ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
    mleFit <- optim(par=log(c(shape.hat,scale.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian=TRUE)
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      scale.hat <- weighted.mean(yi, si)
      s2 <- weighted.mean(yi^2, si) - scale.hat^2
      shape.hat <- scale.hat^2 / s2
      shape.lower <- floor(shape.hat)
      shape.upper <- ceiling(shape.hat)
      shape.hat <- optimize(logProfileLik,
                            lower = shape.lower,
                            upper = shape.upper,
                            maximum = TRUE)$maximum
    } else {
      scale.hat <- weighted.mean(yi, ni)
      s2 <- weighted.mean(yi^2, ni) - scale.hat^2
      shape.hat <- scale.hat^2 / s2
    }
    ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
    mleFit <- optim(par=log(c(shape.hat,scale.hat)),
                    fn=minusLogLik,
                    method="BFGS",
                    hessian=TRUE)
  } ## End of conditional on s.dot == n.dot
  ## estimate <- exp(mleFit$estimate)
  estimate <- exp(mleFit$par)
  newVar <- (1/estimate) %o% (1/estimate)
  observedI <- mleFit$hessian * newVar
  se <- sqrt(diag(solve(observedI)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value


  if (scale) {
    names(estimate) <- c("shape","scale")
    names(se) <- c("shape","scale")
    rFct <- function(shape,scale) -minusLogLik(log(c(shape,scale))) - l
  } else {
    estimate <- c(estimate[1],1/estimate[2])
    names(estimate) <- c("shape","rate")
    se <- se * c(1,estimate[2]^2)
    names(se) <- c("shape","rate")
    rFct <- function(shape,rate) -minusLogLik(log(c(shape,1/rate))) - l
  }

  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)
  
}


llogisMLE <- function(yi,
                      ni = numeric(length(yi))+1,
                      si = numeric(length(yi))+1
                      ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be positive or null")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the location parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
                   )
      cat(txt)
    } else {
      location <- p[1]
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dllogis(yi[si>0],location,scale,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(pllogis(yi[ci>0],location,scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  lyi <- log(yi)
  if (s.dot == n.dot) {
    ## no censored event
    ## Get the empirical mean which is the MLE of parameter mu
    location.moment <- weighted.mean(lyi, si)
    scale.moment <- sqrt(3*(weighted.mean(lyi^2, si) - location.moment^2))/pi
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      location.moment <- weighted.mean(lyi, si)
      scale.moment <- sqrt(3*(weighted.mean(lyi^2, si) - location.moment^2))/pi
    } else {
      location.moment <- weighted.mean(lyi, ni)
      scale.moment <- sqrt(3*(weighted.mean(lyi^2, ni) - location.moment^2))/pi
    }
  }
  ## mleFit <- nlm(minusLogLik,c(location.moment,log(scale.moment)),hessian=TRUE)
  mleFit <- optim(fn=minusLogLik,
                  par=c(location.moment,log(scale.moment)),
                  method="BFGS",
                  hessian=TRUE)
  ## estimate <- c(mleFit$estimate[1],exp(mleFit$estimate[2]))
  estimate <- c(mleFit$par[1],exp(mleFit$par[2]))
  newVar <- c(1,1/estimate[2]) %o% c(1,1/estimate[2])
  se <- sqrt(diag(solve(mleFit$hessian * newVar)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value
  names(estimate) <- c("location","scale")
  names(se) <- c("location","scale")
  rFct <- function(location,scale) -minusLogLik(c(location,log(scale))) - l
  
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}


lnormMLE <- function(yi,
                     ni = numeric(length(yi))+1,
                     si = numeric(length(yi))+1
                     ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be strictly positive")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)
  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the meanlog parameter,\n",
                   "  component 2 is the log of the sdlog parameter.\n"
                   )
      cat(txt)
    } else {
      mean <- p[1]
      sd <- exp(p[2])
      -(ifelse(s.dot>0,sum(dlnorm(yi[si>0],mean,sd,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(plnorm(yi[ci>0],mean,sd,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }
  
  ## Get the empirical mean which is the MLE of parameter mu
  if (s.dot >= 2) {
    mu.hat <- weighted.mean(log(yi), si)
    s2 <- weighted.mean(log(yi)^2, si) - mu.hat^2
  } else {
    mu.hat <- weighted.mean(log(yi), ni)
    s2 <- weighted.mean(log(yi)^2, ni) - mu.hat^2
  }
  if (s.dot == n.dot) {
    ## all events are uncensored
    estimate <- c(mu.hat,sqrt(s2))
    se <- c(sqrt(s2/n.dot),sqrt(s2/(2*n.dot)))
    l <- -minusLogLik(c(estimate[1],log(estimate[2])))
  } else {
    ## numerical fit needed
    ## mleFit <- nlm(minusLogLik,c(mu.hat,0.5*log(s2)),hessian=TRUE)
    mleFit <- optim(fn=minusLogLik,
                    par=c(mu.hat,0.5*log(s2)),
                    method="BFGS",
                    hessian=TRUE)
    ## estimate <- c(mleFit$estimate[1],exp(mleFit$estimate[2]))
    estimate <- c(mleFit$par[1],exp(mleFit$par[2]))
    newVar <- c(1,1/estimate[2]) %o% c(1,1/estimate[2])
    se <- sqrt(diag(solve(mleFit$hessian * newVar)))
    ## l <- -mleFit$minimum
    l <- -mleFit$value
  }
  names(estimate) <- names(formals(dlnorm))[2:3]
  names(se) <- names(formals(dlnorm))[2:3]
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = function(meanlog,sdlog) -minusLogLik(c(meanlog,log(sdlog))) - l,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}

rexpMLE <- function(yi,
                    ni = numeric(length(yi))+1,
                    si = numeric(length(yi))+1
                    ) {

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be positive or null")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the rate parameter,\n",
                   "  component 2 is the log of the rp (refractory period) parameter.\n"
                   )
      cat(txt)
    } else {
      rate <- exp(p[1])
      rp <- exp(p[2])
      -(ifelse(s.dot>0,sum(drexp(yi[si>0],rate=rate,rp=rp,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(prexp(yi[ci>0],rate=rate,
                                    rp=rp,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }
  
  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)
  if (s.dot == 0) stop("No uncensored events")
  rp.hat <- min(yi[si==1]) - .Machine$double.eps
  mu.hat <- weighted.mean(yi-rp.hat, si)
  rate.hat <- 1/mu.hat
  estimate <- c(rate.hat,rp.hat)
  se <- c(rate.hat/sqrt(s.dot),NA)
  l <- -minusLogLik(log(c(rate.hat,rp.hat)))

  names(estimate) <- c("rate","rp")
  names(se) <- c("rate","rp")
  rFct <- function(rate,rp) -minusLogLik(log(c(rate,rp))) - l
  
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}


weibullMLE <- function(yi,
                       ni = numeric(length(yi))+1,
                       si = numeric(length(yi))+1,
                       shape.min = 0.05,
                       shape.max = 5) {
                       
  #####################################################
  ## define internal function makeProfileWeibullLogLik
  #####################################################
  makeProfileWeibullLogLik <- function(data) {

    ## check that data is strictly positive
    if ( any(data <= 0) )
      stop("data should be a vector of positive numbers")
    ## if data is somthing else than a vector (eg a matrix)
    ## coerce it to a vector
    dim(data) <- NULL
    ## get the sample size
    nDot <- length(data)
    ## get the sum of the logs of the variate
    sumLog <- sum(log(data))
    
    logLik <- function(shape) {
      ## check that shape is strictly positive
      if (shape <= 0) {
        return(-Inf)
      } else {
        term1 <- nDot * log(shape)
        term2 <- (shape - 1) * sumLog
        term3 <- - nDot * (1 + log(sum(data^shape)/nDot))
        return( term1 + term2 + term3)
      }
    }

    return(logLik)
  
  }
  ###############################################################
  ## End of internal function makeProfileWeibullLogLik definition
  ###############################################################

  ## check if yi is a spikeTrain object if yes take the "diff"
  if (inherits(yi,"spikeTrain")) yi <- diff(yi)
  ## check that yi elements are positive
  if (any(yi <= 0)) stop("yi elements must be strictly positive")
  ## coerce yi to vector
  yi <- as.numeric(yi)
  ## check that ni has the same length as yi
  if (!identical(length(yi),length(ni))) stop("yi and ni should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(ni < 0)) stop("ni elements must be non-negative")
  if (!identical(ni, round(ni))) stop("ni should be a vector of positive integers")
  if (!identical(length(yi),length(si))) stop("yi and si should have the same length")
  ## check that the elements of ni are non negative and integer
  if (any(si < 0)) stop("si elements must be non-negative")
  if (!identical(si, round(si))) stop("si should be a vector of positive integers")
  if (any(si > ni)) stop("si elements should not be greater than ni elements")

  ## Create function returning the opposite of the log likelihood
  minusLogLik <- function(p) {
    if (missing(p)) {
      txt <- paste("This function argument should be a 2 component vector:\n",
                   "  component 1 is the log of the shape parameter,\n",
                   "  component 2 is the log of the scale parameter.\n"
                   )
      cat(txt)
    } else {
      shape <- exp(p[1])
      scale <- exp(p[2])
      -(ifelse(s.dot>0,sum(dweibull(yi[si>0],shape=shape,scale=scale,log=TRUE)*si[si>0],0))+
        ifelse(c.dot>0,sum(pweibull(yi[ci>0],shape=shape,
                                    scale=scale,lower.tail=FALSE,log.p=TRUE)*ci[ci>0]),0)
        )
    }
  }
  
  ## Get the number of uncensored events
  s.dot <- sum(si)
  ## Get the total number of events
  n.dot <- sum(ni)
  ci <- ni-si
  c.dot <- sum(ci)

  if (s.dot == n.dot) {
    ## no censored event
    mu.hat <- weighted.mean(yi, si)
    s2 <- weighted.mean(yi^2, si) - mu.hat^2
    stat1 <- mu.hat^2/(s2 + mu.hat^2)
    shapeFct <- function(shape) (gamma(1 + 1/shape))^2 / gamma(1 + 2/shape) - stat1
    myRoot <- uniroot(shapeFct, lower = shape.min, upper = shape.max)$root
    shapeInterval <- c(floor(myRoot) - .Machine$double.eps,ceiling(myRoot) + .Machine$double.eps)
    data <- rep(yi[si>0],each=ni[si>0])
    myLogLik <- makeProfileWeibullLogLik(data)
    shape.hat <- optimize(myLogLik, interval = shapeInterval, maximum = TRUE)$maximum
    scale.hat <- (sum(data^shape.hat)/s.dot)^(1/shape.hat)
  } else {
    ## some censored events
    ## if more than 10 events are uncesored get inital guess from them
    ## otherwise use all events
    if (s.dot >= 10) {
      mu.hat <- weighted.mean(yi, si)
      s2 <- weighted.mean(yi^2, si) - mu.hat^2
      stat1 <- mu.hat^2/(s2 + mu.hat^2)
      shapeFct <- function(shape) (gamma(1 + 1/shape))^2 / gamma(1 + 2/shape) - stat1
      myRoot <- uniroot(shapeFct, lower = shape.min, upper = shape.max)$root
      shapeInterval <- c(floor(myRoot) - .Machine$double.eps,ceiling(myRoot) + .Machine$double.eps)
      data <- rep(yi[si>0],each=ni[si>0])
      myLogLik <- makeProfileWeibullLogLik(data)
      shape.hat <- optimize(myLogLik, interval = shapeInterval, maximum = TRUE)$maximum
      scale.hat <- (sum(data^shape.hat)/s.dot)^(1/shape.hat)
    } else {
      mu.hat <- weighted.mean(yi, ni)
      s2 <- weighted.mean(yi^2, ni) - mu.hat^2
      stat1 <- mu.hat^2/(s2 + mu.hat^2)
      shapeFct <- function(shape) (gamma(1 + 1/shape))^2 / gamma(1 + 2/shape) - stat1
      myRoot <- uniroot(shapeFct, lower = shape.min, upper = shape.max)$root
      shapeInterval <- c(floor(myRoot) - .Machine$double.eps,ceiling(myRoot) + .Machine$double.eps)
      data <- rep(yi,each=ni)
      myLogLik <- makeProfileWeibullLogLik(data)
      shape.hat <- optimize(myLogLik, interval = shapeInterval, maximum = TRUE)$maximum
      scale.hat <- (sum(data^shape.hat)/n.dot)^(1/shape.hat)
    }
  } ## End of conditional on s.dot == n.dot
  ## mleFit <- nlm(minusLogLik,log(c(shape.hat,scale.hat)),hessian=TRUE)
  mleFit <- optim(fn=minusLogLik,
                  par=log(c(shape.hat,scale.hat)),
                  method="BFGS",
                  hessian=TRUE)
  ## estimate <- exp(mleFit$estimate)
  estimate <- exp(mleFit$par)
  newVar <- (1/estimate) %o% (1/estimate)
  observedI <- mleFit$hessian * newVar
  se <- sqrt(diag(solve(observedI)))
  ## l <- -mleFit$minimum
  l <- -mleFit$value

  names(estimate) <- c("shape","scale")
  names(se) <- c("shape","scale")
  rFct <- function(shape,scale) -minusLogLik(log(c(shape,scale))) - l
  
  result <- list(estimate = estimate,
                 se = se,
                 logLik = l,
                 r = rFct,
                 mll = minusLogLik,
                 call = match.call()
                 )
  class(result) <- "durationFit"
  return(result)

}

is.durationFit <- function(obj) {

  if (!("durationFit" %in% class(obj))) return(FALSE)
  
  expectedNames <- c("estimate",
                     "se",
                     "logLik",
                     "r",
                     "mll",
                     "call")
  if (!all(expectedNames %in% names(obj))) return(FALSE)
  TRUE
    
}

qqDuration <- function(durationFit,
                       CI=c(0.95,0.99),
                       type="l",
                       xlab,
                       ylab,
                       main,
                       sub,
                       ylim,
                       dataLwd=2,
                       ablineCol=2,
                       ...
                       ) {
  ## Check that durationFit is a durationFit object
  if (!is.durationFit(durationFit))
    stop(paste(deparse(substitute(durationFit)),
               "should be a durationFit object."
               )
         )

  ## check CI
  if (!is.null(CI)) {
    ## make sure that CI is at most of length 2 otherwise
    ## keep the first 2 components
    if (length(CI) > 2) CI <- CI[1:2]
    ## Check that each component of CI is between 0 and 1
    if (any(CI>=1 | CI<=0))
      stop(paste(deparse(substitute(CI)),
                 "components should be in (0,1)")
           )
  } ## End of conditional on !is.null(CI)
  
  ## get the uncensored events
  ## The next code line causes an unimportant warning during
  ## automatic code check by "R CMD check"
  x <- evalq(rep(yi[si>0],each=si[si>0]),envir=environment(durationFit$r))
  ## get the y values of the QQ plot
  y <- sort(x)
  ## get the theoretical x values
  distribution <- switch(deparse(durationFit$call[[1]]),
                         invgaussMLE = qinvgauss,
                         gammaMLE = qgamma,
                         weibullMLE = qweibull,
                         lnormMLE = qlnorm,
                         rexpMLE = qrexp,
                         llogisMLE = qllogis)

  if (is.null(distribution))
    stop(paste("Unknown distribution in",
               deparse(substitute(durationFit))
               )
         )

  sampleSize <- length(y)
  pp <- ppoints(sampleSize)
  x.theo <- do.call(distribution,
                    c(list(p=pp),
                      as.list(durationFit$estimate)
                      )
                    )


  if (!is.null(CI)) {
    ## Get the confidence intervals
    ci <- sapply(1:sampleSize,
                 function(idx)
                 as.vector(sapply(CI,
                                  function(l) {
                                    thePs <- c((1-l)/2,1-(1-l)/2)
                                    thePs <- qbeta(thePs,
                                                   idx,
                                                   sampleSize-idx+1)
                                    do.call(distribution,c(list(p=thePs),
                                                           as.list(durationFit$estimate))
                                            )
                                    
                                  }
                                  )
                           )
                 )
  } ## End of conditional on !is.null(CI)
  
  modelName <- switch(deparse(durationFit$call[[1]]),
                      invgaussMLE = "invgauss",
                      gammaMLE = "gamma",
                      weibullMLE = "weibull",
                      lnormMLE = "lnorm",
                      rexpMLE = "rexp",
                      llogisMLE = "llogis")
  sampleName <- deparse(durationFit$call[["yi"]])

  if (missing(xlab)) xlab <- paste("Quantiles of",
                                   modelName)
  if (missing(ylab)) ylab <- "Sample quantiles"
  if (missing(main)) main <- paste("QQ plot of",
                                   sampleName,
                                   "vs fitted",
                                   modelName)
  if (missing(sub)) {
    sub <- paste(sampleSize,
                 "uncensored intervals used.")
    if (!is.null(CI))
      sub <- paste(sub,
                   "CI at:",
                   paste(CI,collapse=",")
                   )
  }

  if (missing(ylim)) ylim <- c(min(min(x.theo),min(y)),
                               max(max(x.theo),max(y))
                               )
  
  plot(x.theo,y,
       xlab=xlab,
       ylab=ylab,
       main=main,
       sub=sub,
       ylim=ylim,
       type="n",
       ...
       )
  abline(a=0,b=1,col=ablineCol)
  if (!is.null(CI)) apply(ci,1,function(x) lines(x.theo,x,lty=2))
  lines(x.theo,y,type=type,lwd=dataLwd)
  
}

compModels <- function(yi,
                       ni = numeric(length(yi))+1,
                       si = numeric(length(yi))+1,
                       models = c("invgauss","lnorm","gamma","weibull","llogis","rexp"),
                       type = c("qq","s"),
                       log = TRUE,
                       plot=TRUE
                       ) {

  if (is.spikeTrain(yi)) yi <- diff(yi)
  ## Check that the models are known
  knownModels <- c("invgauss",
                   "lnorm",
                   "gamma",
                   "weibull",
                   "llogis",
                   "rexp")

  if (!all(models %in% knownModels))
    stop("Some requested models are not implemented.")

  theFits <- lapply(models,
                    function(m) switch(m,
                                       invgauss = invgaussMLE(yi,ni,si),
                                       lnorm = lnormMLE(yi,ni,si),
                                       gamma = gammaMLE(yi,ni,si),
                                       weibull = weibullMLE(yi,ni,si),
                                       llogis = llogisMLE(yi,ni,si),
                                       rexp = rexpMLE(yi,ni,si)
                                       )
                    )

  if (plot) {
    myNrow <- (length(models)-1) %/% 2 + 1
    layout(matrix(1:(myNrow*2),nrow=myNrow))
    oldpar <- par(mar=c(4,3,2,1))
    on.exit(par(oldpar))
    if (type[1] == "qq") {
      sapply(1:length(models),
             function(idx) {
               if (log) {
                 qqDuration(theFits[[idx]],
                          main="",
                            ylab="",
                            sub="",
                            log="xy")
               } else {
                 qqDuration(theFits[[idx]],
                            main="",
                            ylab="",
                            sub="")
               }
             }
             )
      
    } else {
      
      isi <- rep(yi[ni>0],each=ni[ni>0])
      status <- unlist(lapply(1:length(ni),
                              function(idx) {
                                if (ni[idx]==0) return(numeric())
                                else {
                                  result <- numeric(ni[idx])
                                  if (si[idx]>0) result[1:si[idx]] <- 1
                                  return(result)
                                }
                              }
                              )
                       )
      kmFit <- survfit(Surv(isi,status) ~1)
      sapply(1:length(models),
             function(idx) {
               plot(kmFit,log=log,mark.time=FALSE,main=models[idx])
               survFct <- switch(models[idx],
                                 invgauss = pinvgauss,
                                 lnorm = plnorm,
                                 gamma = pgamma,
                                 weibull = pweibull,
                                 llogis = pllogis,
                                 rexp = prexp)
               isi <- sort(isi)
               y <- do.call(survFct,c(list(q=isi,lower.tail=FALSE),as.list(theFits[[idx]]$estimate)))
               lines(isi,y,col=2)
             }
             )
    }
  } ## End of conditional on plot
  aic <- sapply(theFits, function(f) -2*f$logLik+4)
  ordered <- sort.int(aic,index.return=TRUE)
  result <- ordered$x
  names(result) <- models[ordered$ix]
  result

}

dinvgauss <- function(x,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      log = FALSE
                      ){
  
  if( any(x <= 0) ) stop("y must contain positive values")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  tmp <- -(x-mu)^2/(2*x*sigma2*mu^2)-(log(2*pi*sigma2)+3*log(x))/2
  if(!log) tmp <- exp(tmp)
  tmp

}

dllogis <- function(x,
                    location = 0,
                    scale = 1,
                    log = FALSE) {
  
  ## Check that x elements are strictly positive
  if (any(x <= 0)) stop("x elements must be strictly positive")

  if (!log) dlogis(log(x),location,scale)/x
  else dlogis(log(x),location,scale,log=TRUE) - log(x)
  
}

drexp <- function(x,
                  rate = 10,
                  rp = 0.005,
                  log = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  x <- x-rp
  dexp(x,rate=rate,log=log)
}


hgamma <- function(x,
                   shape,
                   rate = 1,
                   scale = 1/rate,
                   log = FALSE) {

  if (!log) dgamma(x,shape,,scale)/pgamma(x,shape,,scale,lower.tail=FALSE)
  else dgamma(x,shape,,scale,log=TRUE) -
    pgamma(x,shape,,scale,lower.tail=FALSE,log.p=TRUE)

}

hinvgauss <- function(x,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      log = FALSE) {

  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  
  t <- x/mu
  v <- sqrt(x*sigma2)
  cutOff <- 1e-12
  bigEnough <- pinvgauss(x, mu, sigma2, lower.tail = FALSE) > cutOff
  logIntensity <- -( (t[bigEnough]-1)^2/(x[bigEnough]*sigma2) + log(2*sigma2*pi*x[bigEnough]^3) )/2 -
    log(1-pnorm((t[bigEnough]-1)/v[bigEnough])-exp(2/(mu*sigma2))*pnorm(-(t[bigEnough]+1)/v[bigEnough])
        )
  if (sum(bigEnough) < length(x)) {
    cutOffQ <- qinvgauss(1 - cutOff, mu, sigma2)
    cutOffValue <- dinvgauss(cutOffQ, mu, sigma2, log = TRUE) - log(cutOff)
    logIntensity <- c(logIntensity,rep(cutOffValue,length(x) - sum(bigEnough))) 
  }
  if (!log) return(exp(logIntensity))
  else return(logIntensity)
  
}

hllogis <- function(x,
                    location = 0,
                    scale = 1,
                    log = FALSE) {
  
  ## Check that x elements are strictly positive
  if (any(x <= 0)) stop("x elements must be strictly positive")

  if (!log) pllogis(x,location,scale)/(scale*x)
  else pllogis(x,location,scale,log.p=TRUE) - log(scale) - log(x)

}

hlnorm <- function(x,
                   meanlog = 0,
                   sdlog = 1,
                   log = FALSE) {
  ## Check that y elements are strictly positive
  if (any(x <= 0))
    stop(paste("The elements of",
               deparse(substitute(x)),
               "should be strictly positive.")
         )

  if (!log) dlnorm(x,meanlog,sdlog)/plnorm(x,meanlog,sdlog,lower.tail=FALSE)
  else dlnorm(x,meanlog,sdlog,log=TRUE)-
    plnorm(x,meanlog,sdlog,lower.tail=FALSE,log.p=TRUE)
  
}

hrexp <- function(x,
                  rate = 10,
                  rp = 0.005,
                  log = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  
  if (!log) ifelse(x >= rp, rate, 0)
  else ifelse(x >= rp, log(rate), -Inf)
}


hweibull <- function(x,
                     shape,
                     scale = 1,
                     log = FALSE) {

  if (!log) dweibull(x, shape, scale) / pweibull(x,shape,scale,lower.tail=FALSE)
  else dweibull(x,shape,scale,log=TRUE) - pweibull(x,shape,scale,lower.tail=FALSE,log.p=TRUE)
  
}

isiHistFit <- function(spikeTrain,
                       model,
                       nbins = 10,
                       CI=0.95,
                       ...) {

  spikeTrainName <- deparse(substitute(spikeTrain))
  if (!is.spikeTrain(spikeTrain)) spikeTrain <- as.spikeTrain(spikeTrain)

  isi <- diff(spikeTrain)
  knownModels <- c("invgauss", "lnorm", "gamma", "weibull", "llogis", "rexp")
  model <- model[1]
  if (!model %in% knownModels) 
    stop(paste(deparse(substitute(model)),
               "is not implemented.")
         )

  CI <- CI[1]
  if (CI <= 0 | CI >= 1)
    stop(paste("A CI of", CI, "does not make sense."))

  theFit <- switch(model,
                   invgauss = invgaussMLE(isi),
                   lnorm = lnormMLE(isi),
                   gamma = gammaMLE(isi),
                   weibull = weibullMLE(isi),
                   llogis = llogisMLE(isi),
                   rexp = rexpMLE(isi)
                   )

  pSeq <- (1:(nbins-1))/nbins

  theEstimates <- theFit$estimate
  qFct <- switch(model,
                 invgauss = function(p) qinvgauss(p,theEstimates[1],theEstimates[2]),
                 lnorm = function(p) qlnorm(p,theEstimates[1],theEstimates[2]),
                 gamma = function(p) qgamma(p,theEstimates[1],theEstimates[2]),
                 weibull = function(p) qweibull(p,theEstimates[1],theEstimates[2]),
                 llogis = function(p) qllogis(p,theEstimates[1],theEstimates[2]),
                 rexp = function(p) qrexp(p,theEstimates[1],theEstimates[2])
                 )

  breaks <- c(0,qFct(pSeq),max(isi)+.Machine$double.eps)
  mids <- breaks[-(nbins+1)] + diff(breaks)/2
  main <- paste("Isi histogram and fitted",
                model,
                "distribution for",
                spikeTrainName
                )
  sub <- paste("CI at ",
               CI*100,
               "%. Sample size: ",
               length(isi),".",
               sep=""
               )
  hist(isi,
       breaks=breaks,
       main=main,
       sub=sub,
       xlab="isi (s)",
       ...)

  X <- seq(.Machine$double.eps,breaks[nbins+1],length.out=501)

  dFct <- switch(model,
                 invgauss = function(x) dinvgauss(x,theEstimates[1],theEstimates[2]),
                 lnorm = function(x) dlnorm(x,theEstimates[1],theEstimates[2]),
                 gamma = function(x) dgamma(x,theEstimates[1],theEstimates[2]),
                 weibull = function(x) dweibull(x,theEstimates[1],theEstimates[2]),
                 llogis = function(x) dllogis(x,theEstimates[1],theEstimates[2]),
                 rexp = function(x) drexp(x,theEstimates[1],theEstimates[2])
                 )
  Y <- dFct(X)
  lines(X,Y,col=2)
  ci <- qbinom(c((1-CI)/2,1-(1-CI)/2),size=length(isi),prob=1/nbins)
  l <- ci[1]/(length(isi)*diff(breaks))
  u <- ci[2]/(length(isi)*diff(breaks))
  invisible(sapply(1:nbins,function(idx) segments(mids[idx],l[idx],mids[idx],u[idx],col=2)))
  
}

pinvgauss <- function(q,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL,
                      lower.tail = TRUE,
                      log.p = FALSE
                      ){
  
  if( any(q <= 0) ) stop("q must contain positive values")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  t <- q/mu
  v <- sqrt(q*sigma2)

  ## Use Eq. 4 of Whitemore GA and Seshadri V (1987)
  ## The American Statistician 41:280-281
  if (lower.tail & !log.p)
    return(pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v))
  if (!lower.tail & !log.p)
    return(1 - (pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v)))
  if (lower.tail & log.p)
    return(log(pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v)))
  if (!lower.tail & log.p)
    return(log(1 - (pnorm((t-1)/v)+exp(2/(mu*sigma2))*pnorm(-(t+1)/v))))
  
}

pllogis <- function(q,
                    location = 0,
                    scale = 1,
                    lower.tail = TRUE,
                    log.p = FALSE) {
  
  plogis(log(q),location,scale,lower.tail,log.p)
  
}

prexp <- function(q,
                  rate = 10,
                  rp = 0.005,
                  lower.tail = TRUE,
                  log.p = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  q <- q-rp
  pexp(q,rate,lower.tail,log.p)
}

qinvgauss <- function(p,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL
                      ){

  if( any(p < 0 | p > 1) ) stop("p must lie between 0 and 1")
  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }
  
  len <- max(length(p),length(mu),length(sigma2))

  if(length(p) != len) {
    if(length(p) == 1) p <- rep(p,len)
    else stop("length of p incorrect")
  }
  if(length(mu) != len) {
    if(length(mu) == 1) mu <- rep(mu,len)
    else stop("length of m incorrect")
  }
  if(length(sigma2) != len) {
    if(length(sigma2) == 1) sigma2 <- rep(sigma2,len)
    else stop("length of sigma2 incorrect")
  }

  ## Use Whitemore and Yalovky (1978, Technometrics, 20:207-208)
  ## approximation to get starting value for the numerical
  ## inversion of the cumulative distribution function.
  theta <- 1/mu/sigma2
  approx <- mu * exp(qnorm(p)*sqrt(1/theta)-0.5/theta)
  sapply(1:len, function(idx) {
    if (identical(p[idx],0)) return(0)
    if (identical(p[idx],1)) return(Inf)
    interval <- approx[idx]*c(0.95,1.05)
    h <- function(q) pinvgauss(q, mu[idx], sigma2[idx]) - p[idx]
    while (h(interval[1])*h(interval[2]) > 0)
      interval <- interval*c(0.9,1.1)
    uniroot(h,interval)$root
  }
         )
  
}

qllogis <- function(p,
                    location = 0,
                    scale = 1,
                    lower.tail = TRUE,
                    log.p = FALSE) {
  
  exp(qlogis(p,location,scale,lower.tail,log.p))
  
}

qrexp <- function(p,
                  rate = 10,
                  rp = 0.005,
                  lower.tail = TRUE,
                  log.p = FALSE) {

  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  rp + qexp(p,rate,lower.tail,log.p)

}


rinvgauss <- function(n = 1,
                      mu = 1,
                      sigma2 = 1,
                      boundary = NULL
                      ){

  if( any(mu <= 0) ) stop("mu must be positive")
  if (all(!is.null(sigma2)))
    if( any(sigma2 <= 0) ) stop("sigma2 must be positive")
  if (all(!is.null(boundary)))
    if( any(boundary <= 0) ) stop("boundary must be positive")
  if ( all(!is.null(sigma2)) && all(!is.null(boundary)) )
    stop("One of sigma2 or boundary must be specified, not both")

  if (all(!is.null(boundary))) {
    ## We work with the parameterization in term of boundary (ie, sigma2 = 1)
    ## We convert it in term of mu and sigma2
    sigma2 <- (1/boundary)^2
    mu <- boundary*mu
  }

  ## Use method of Michael JR, Schucany WR and Haas RW (1976)
  ## The American Statistician, 30:88-90
  v0 <- rchisq(n,1)
  x1 <- mu + 0.5*mu^2*sigma2*v0 - 0.5*mu*sigma2*sqrt(4*mu*v0/sigma2+mu^2*v0^2)
  ifelse(rbinom(length(x1),1,mu/(mu+x1)) == 1,x1,mu^2/x1)
}


rllogis <- function(n,
                    location = 0,
                    scale = 1) {
  
  exp(rlogis(n,location,scale))

}

rrexp <- function(n,
                  rate = 10,
                  rp = 0.005) {
  ## Check that rate is strictly positive
  if (rate <= 0) stop("rate should be strictly positive")
  ## Check that rp is positive or null
  if (rp < 0) stop("rp should be positive or null")
  rp+rexp(n,rate)
}

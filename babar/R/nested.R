########################################################################
# nested.R - code to estimate posterior and evidence via nested sampling
# (c) Matthew Hartley, Nick Pullen and Richard Morris
########################################################################
##source("params2priors.R")

.logPlus <- function(a, b) {
  # Sum a and b via their logarithms, such that if:
  # a = log(x) and b = log(y)
  # logPlus(a, b) returns log(x + y)
  if (a > b) {
    sum = a + log(1+exp(b-a))
  } else {
    sum = b + log(1+exp(a-b))
  }
  return(sum)
}

.makeStep <- function(current.values, step.vector) { # Make MC step
  # Make random step in parameter space
  #
  # Args:
  #   current.values: Vector of current parameter values.
  #   step.vector: Vector of step sizes, of the same length as the parameter
  #                values.
  #
  # Returns:
  #   New values as a vector of length number.of.parameters
  n <- length(step.vector)
  
  stopifnot(length(current.values) == length(step.vector))
  
  ## print(c(current.values, step.vector, n, runif(n, -1, 1)))
  return(current.values + step.vector * runif(n, -1, 1))
}

.makeBoundedStep <- function(current.values, step) {
  # Make step in parameter space. If step exceeds bounds, instead pick point in
  # space randomly chosen using uniform distribution.
  #
  # Args:
  #   current.values: Vector of parameter values
  #   
  step.attempt <- .makeStep(current.values, step)
  
  ## print(c(current.values, step, step.attempt))
  stopifnot(length(current.values) == length(step))## Note to self: Need to make length(lower.bounds) depend on numParams
  ##  failed.steps <- step.attempt < lower.bounds | step.attempt > upper.bounds
  failed.steps <- step.attempt < 0 | step.attempt > 1## Bounds of unit hypercube
  
  for (i in which(failed.steps)) {
    step.attempt[i] <- runif(1, min = 0, max = 1)
    ## DEBUG     cat("***HERE*******\n")
  }
  
  # Replace value outside the bounds by randomly chosen point in parameter space
  #  if((step.attempt < lower.bounds) + (step.attempt > upper.bounds) > 0) {
  #   step.attempt = runif(1, lower.bounds, upper.bounds)
  #}
  #step.attempt[which(outside.bounds)] = runif(1, lower.bounds, upper.bounds)
  return(step.attempt)
}

.makeBoundedSingleStep <- function(current.value, step) {
  # Make a step of single (scalar) parameter value. If step takes us outside bounds,
  # pick point in space from uniform distribution within bounds
  #
  # Args:
  #   current.value: scalar current parameter value
  #   step: scalar step size
  #
  # Returns:
  #   New parameter value
  
  step.attempt <- current.value + step * runif(1, -1, 1)
  
  if (step.attempt < 0 | step.attempt > 1) {
    step.attempt <- runif(1, min = 0, max = 1)
  }
  
  return(step.attempt)
}

.ballMcmcStep <- function( # Make MCMC step in parameter space for all parameters simultaneously
  ### Make single MCMC step from current.values as follows:
  ### 1) Attempt step in all parameter directions simultaneously
  ### 2) If step is successful (ll bigger than llMin), note this
  current.values, 
  ### Vector of current parameter values.
  llFun, 
  ### Function to calculate loglikelihood from llFun(param.vector).
  llMin, 
  ### Current minimum loglikelihood.
  steps 
  ### Vector of step values.
) {
  candidate <- .makeBoundedStep(current.values, steps)
  
  ## print(c(candidate, "********", current.values, steps))
  if (llFun(candidate) > llMin) {
    accepted <- TRUE
    new.values <- candidate
  } else {
    accepted <- FALSE
    new.values <- current.values
  }
  
  return(list(new.values=new.values, accepted=accepted))
  ###    accepted: a vector of boolean values describing whether each
  ###    parameter's step was accepted
  ###
  ###    new.values: a vector of the new parameter values after the step
}

.makeMcmcStep <- function( # Make MCMC step in parameter space one parameter at a time
  ### Make single MCMC step from current.values as follows:
  ### 1) Attempt step in each parameter direction
  ### 2) If step is successful (ll bigger than llMin), note this
  current.values, 
  ### current.values: Vector of current parameter values.
  llFun, 
  ### llFun: Function to calculate loglikelihood from
  ### llFun(param.vector).
  llMin,
  ### llMin: Current minimum loglikelihood.
  steps                     
  ### steps: Vector of step values.
) { 
  n <- length(current.values)
  accepted <- rep(FALSE, n)
  for (i in 1:n) {
    candidate <- current.values
    # Attempt step along one axis in parameter space
    candidate[i] <- .makeBoundedSingleStep(current.values[i], steps[i])
    if(llFun(candidate) > llMin) {
      current.values <- candidate
      accepted[i] = TRUE
    }
  }
  return(list(accepted=accepted, new.values=current.values))
  ###    accepted: a vector of boolean values describing whether each
  ###    parameter's step was accepted
  ###
  ###    new.values: a vector of the new parameter values after the step
}

ballExplore <- structure(function
### Explore around the current values in a sphere in parameter space
(current.values,
### Vector of current parameter values
 steps,
### Vector of steps for current parameters
 llMin,
 llFun
 ) {
  
  msub <- 20
  m <- 5
  
  ## print(c(current.values, steps, llMin, llFun))
  accepted <- vector()
  for (k in 1:m) {
    for (i in 1:msub) {
      ret <- .ballMcmcStep(current.values, llFun, llMin, steps)
      current.values <- ret$new.values
      accepted <- cbind(accepted, ret$accepted)
    }
  }
  
  ## print(c(ret$new.values, llFun(ret$new.values), "^^^^^^^^^^^^^", accepted))
  ratio <- sum(accepted) / length(accepted)
  
  if (ratio > 0.6) steps <- steps * (1 + 2 * (ratio - 0.6) / 0.4)
  if (ratio < 0.4) steps <- steps / (1 + 2 * ((0.4 - ratio) / 0.4))
  
  steps <- pmax(steps, 0.01)## Feel this should perhaps alyways just depend on steps
  
  return(list(new.values=ret$new.values, new.steps=steps, new.ll=llFun(ret$new.values)))
})

.explore <- function(current.values, steps, llMin, llFun) {
  # Explore parameter space around supplied point
  # returns new point and new step size
  msub <- 20
  m <- 5
  
  for(k in 1:m) {
    accepted <- vector()
    for(i in 1:msub) {
      ret <- .makeMcmcStep(current.values, llFun, llMin, steps)
      current.values <- ret$new.values
      accepted <- cbind(accepted, ret$accepted)
    }
    
    # DEBUG
    #cat(sum(accepted[1,]), " ", sum(accepted[2,]), "\n")
    
    for(i in 1:length(current.values)) {
      ratio <- sum(accepted[i,]) / length(accepted[i,])
      if (ratio > 0.6) {
        steps[i] <- steps[i] * (1 + 2 * (ratio - 0.6) / 0.4)
      }
      if (ratio < 0.4) {
        steps[i] <- steps[i] / (1 + 2 * ((0.4 - ratio) / 0.4))
      }
      steps[i] <- max(steps[i], 0.01) ## Feel this should perhaps alyways just depend on steps
      
      #cat(sprintf("%f, %f\n", ratio, step[i]))
    }
  }
  return(list(new.values=current.values, new.step=steps, new.ll=llFun(current.values)))
}

.generatePriorSamples <- function( # Generate prior samples
  ### Generate a set of prior samples from a unit hypercube as recommended by Skilling
  n.samples, 
  ### The number of prior samples to be generated.
  numberOfParameters
  ### The number of parameters to be inferred
)  {
  prior.samples <- replicate(n.samples,  .generateUnitHyperCube(numberOfParameters))
  
  ## If we only have a 1d problem, we need to ensure that we return a matrix of
  ## suitable dimensions, by default we end up with a vector, where ncol returns
  ## NULL, so we use t() to turn it into a suitably dimensioned matrix
  if (numberOfParameters == 1) {
    return(as.matrix(prior.samples))
  } else {
    return(t(prior.samples))
  }
  ###   Matrix of prior samples, of dimension n.samples x n.params.
}

.generateUnitHyperCube <- function( # Generate U(0,1) samples
  ### Generates n samples drawn from the uniform distribution U(0,1)
  n
  ### The number of parameters in the problem
) {
  prior <- runif(n, min = 0, max = 1)
  
  return(prior)
  ###   A vector of values drawn from the uniform distribution
}

nestedSampling <- structure(function
  ### Perform nested sampling for the given model (expressed through the log
  ### likelihood function).
  (llFun,
  ### The loglikelihood function, which should be of the form:
  ### llFun(params) where params is a vector containing a single
  ### instance of parameters, and should return the log likelihood for
  ### those parameters.
  numberOfParameters,
  ### The number of parameters to be inferred. Should always equal the number of parameters in the likelihood function.
  prior.size,
  ### The number of initial prior objects to sample. This is the size of the active set of samples maintained throughout the procedure.
  transformParams,
  ### Function to transform parameters from the unit hypercube into the correct distribution
  exploreFn=ballExplore,
  ### Preferred exploration routine. Defaults to a dodgy implementation
  ### currently but ballExplore likely to be more robust.
  tolerance=0.1
  ### The tolerance at which the routine should stop. Changing this will
  ### affect how long it takes to run nested sampling. Larger values
  ### will terminate sooner but may give inaccurate evidence
  ### scores. This is the expected maximum remaining contribution to
  ### logZ from the points in the active set. delta Z_i = L_max*X_i <
  ### tolerance.
) {
  
  prior.samples <- .generatePriorSamples(prior.size, numberOfParameters)
  # Calculate the log likelihoods for the prior samples
  ll.values <- apply(prior.samples, 1, llFun)
  evaluated.samples = cbind(ll.values, prior.samples)
  ordered.samples <- evaluated.samples[order(evaluated.samples[,1]),]
  
  numParamsAddOne = NCOL(ordered.samples)
  if (numParamsAddOne == 1) { numParamsAddOne = 2; ordered.samples = t(ordered.samples)} ## Stupid hack to get round the fact that R can't ever treat a vector as having 2 columns and 1 row. Quite an edge case because only rears its head when 1 prior object is chosen.
  numParams <- numParamsAddOne - 1L
  posterior.samples <- matrix(ncol = numParamsAddOne)
  first <- TRUE
  
  steps <- rep(1.0, numParams) 
  
  # Initialise parameters for evidence calculation
  logZ <- -2.2e308
  log.width <- log(1 - exp(-1 / prior.size))
  log.weight = numeric()## This is going to be inefficient later
  entropy <- 0
  maxLeft = Inf 
  logTolerance = log(tolerance)
  samplesCounter = 0
  
  while(maxLeft > logTolerance) {
    samplesCounter = samplesCounter + 1
    ## Add the worst sample to the posterior samples
    if ( first ) {
      ##      worst <- ordered.samples[1,]
      posterior.samples[1,] <- ordered.samples[1,]
      first <- FALSE
    } else {
      posterior.samples <- rbind(posterior.samples, ordered.samples[1,])
    }
    
    ## Update evidence
    ## llMin <- min(ordered.samples[,1])
    llMin <- ordered.samples[1,1]
    log.weight[samplesCounter] = llMin + log.width
    logZnew <- .logPlus(logZ, log.weight[samplesCounter])
    ## worst.weight <- log.width + worst[1]
    worst.weight <- log.weight[samplesCounter]
    ## print(c("*********", llMin, worst.weight, logZnew))
    ## print(entropy)
    if (is.infinite(logZ)) {
      entropy <- exp(worst.weight - logZnew) * llMin - logZnew ## Hopefully this is what Skilling intended if your language doesn't allow 0*Inf==0
    } else {
      entropy <- exp(worst.weight - logZnew) * llMin + exp(logZ - logZnew) * (entropy + logZ) - logZnew
    }
    #stop()
    logZ <- logZnew
    log.width = log.width - 1 / prior.size
    # DEBUG cat(sprintf("llMin: %f, ev: %f\n", llMin, logZ))
    ## print(c(llMin, logZnew, entropy, worst.weight, logZ, log.width))
    ## Randomly select a sample that isn't the worst
    if (numParams == 1) { ## Hope people don't choose a prior of size 1. As this is another thing that needed changing.
      selected.point <- 1
    } else {
      selected.point <- sample(2:prior.size, 1)
    }    
    llSelected <- ordered.samples[selected.point,1]
    
    ## Explore around that point, updating step size as we go
    ret <- exploreFn(ordered.samples[selected.point, 2:numParamsAddOne], 
                     steps, llMin, llFun)
    steps <- ret$new.step
    new.point <- ret$new.values
    new.ll <- ret$new.ll
    ## print(c(selected.point, llSelected, steps, new.point, new.ll))
    #stop()
    ## Replace the worst sample by the results of explore()
    ordered.samples[1,] <- c(new.ll, new.point)
    
    ## Re-sort the samples (not very efficient, but we don't care)
    ordered.samples <- ordered.samples[order(ordered.samples[,1]),]
    ## Update max loglikelihood and remaining prior volume to check termination
    if (numParams == 1 & prior.size == 1) { ordered.samples = t(ordered.samples) } ## Hope people don't choose a prior of size 1. As this is another thing that needed changing.
    max.ll <- ordered.samples[prior.size, 1]
    remaining.prior.volume <- exp(-samplesCounter / prior.size)
    deltaZi = max.ll + log(remaining.prior.volume)
    maxLeft = deltaZi - logZ
    #DEBUG cat(sprintf("maxll: %f, logRemPriorVol: %f, their sum: %f, logZrem: %f, sampleNum: %d \n", max.ll, log(remaining.prior.volume), deltaZi, deltaZi - logZ, samplesCounter))
    
    ## print(c(entropy, sqrt(entropy)))
    # DEBUG print(ordered.samples[2:3,prior.size])
  }
  ## Now we add the final contribution to logZ from the active set
  # DEBUG cat("logZ before final correction",logZ,"\n")
  lw <- log(remaining.prior.volume/prior.size)
  for (i in 1:prior.size) {
    lL <- ordered.samples[i, 1]
    log.weight[i + samplesCounter] = lw + lL
    logZ <- .logPlus(logZ, log.weight[i + samplesCounter])
  }
  
  logZerror <- sqrt(entropy / prior.size)
  # DEBUG cat("logZerror",logZerror,"entropy",entropy,"\n")
  # DEBUG cat("logweight",log.weight,"len",length(log.weight),"\n")
  
  posterior.samples = rbind(posterior.samples, ordered.samples)
  # DEBUG cat(nrow(posterior.samples),"****", ncol(posterior.samples),"\n")
  
  ## Transform the unit parameters to their true values
  if (numParams==1) {
    print("This is experimental if numParams = 1. Needs more testing!")
    posterior.samples = cbind(posterior.samples[,1], sapply(posterior.samples[,-1], transformParams))
  } else {
    ## For some reason apply returns the transpose  of what I need. Grrr!
    posterior.samples[,-1] = t(apply(posterior.samples[,-1], 1, transformParams))
  }
  
  ## posterior.samples now becomes numParams + 2 columns, first is logWeights, 2nd is logLikes then the true parameter values
  posterior.samples = cbind(log.weight, posterior.samples)
  
  ## Calculating the means and variances over each parameter
  totalSamples = nrow(posterior.samples)
  sampleMean <- rep(0,numParams)
  sampleMeanOfSquare <- rep(0,numParams)
  sampleVariance <- rep(0,numParams)
  
  for (i in 1:totalSamples) {
    sampleWeight = exp(posterior.samples[i,1] - logZ)## normalising
    for (j in 1:numParams) {
      sampleMean[j] = sampleMean[j] + sampleWeight*posterior.samples[i, j+2]
      sampleMeanOfSquare[j] = sampleMeanOfSquare[j] + sampleWeight*posterior.samples[i, j+2]*posterior.samples[i, j+2]
    }
  }
  for (k in 1:numParams) {
    if (sampleMeanOfSquare[k] < sampleMean[k]*sampleMean[k]) {
      sampleMeanOfSquare[k] = 0.0
    } else {
      sampleVariance[k] = sampleMeanOfSquare[k] - sampleMean[k]*sampleMean[k] 
    }
  }
  
  colnames(posterior.samples) = c("logWeight", "logLikelihood", paste0("parameter",seq(2, numParamsAddOne) -1))
  
  return(list(logevidence=unname(logZ), posterior=posterior.samples, logZerror = unname(logZerror), parameterMeans = sampleMean, parameterVariances = sampleVariance, entropy = unname(entropy)))
  ### logevidence: The logarithm of the evidence, a scalar.
  ###
  ### posterior.samples: The samples from the posterior, together with their log weights and 
  ###                      log likelihoods as a m x n matrix, where m is the
  ###                      number of posterior samples and n is the number of
  ###                      parameters + 2. The log weights are the first column and the log likelihood values are the
  ###                      second column of this matrix. The sum of the log-weights = logZ.
  ###
  ### logZerror: An estimate of the numerical uncertainty of the log evidence.
  ###
  ### parameterMeans: A vector of the mean of each parameter, length = no. of parameters.
  ###
  ### parameterVariances: A vector of the variance of each parameter, length = no. of parameters.
  ###
  ### entropy: The information --- the natural logarithmic measure of the prior-to-posterior shrinkage.
}, ex=function() {

  mu <- 1
  sigma <- 1
  data <- rnorm(100, mu, sigma)
  
  transform <- function(params) {
    tParams = numeric(length=length(params))
    tParams[1] = GaussianPrior(params[1], mu, sigma)
    tParams[2] = UniformPrior(params[2], 0, 2 * sigma)
    return(tParams)
  }

  llf <- function(params) {
    tParams = transform(params)
    mean = tParams[1]
    sigma = tParams[2]
    n <- length(data)
    ll <- -(n/2) * log(2*pi) - (n/2) * log(sigma**2) - (1/(2*sigma**2)) * sum((data-mean)**2)
    return(ll)
  }
  
  prior.size <- 25
  tol <- 0.5

  ns.results <- nestedSampling(llf, 2, prior.size, transform, tolerance=tol)

})

.staircaseSampling = function(
  ### Generate a set of equally weighted posterior samples from the output of nested sampling  
  posterior
  ### Matrix from the output of nested sampling. Should have log weights
  ### in column 1, log likelihoods in column 2 and then the parameter
  ### values in the remaining column(s).
) {
  logSamples = 0.0
  logevidence = log(sum(exp(posterior[,1])))## check then remove ret$'s below
  numNestedSamples = nrow(posterior)
  
  for (i in 1:numNestedSamples) { ## look to vectorise using sum perhaps
    logSamples = logSamples -
      exp(posterior[i,1] - logevidence)*(posterior[i,1] - logevidence)
  }
  
  multiplicity = exp(posterior[,1] - logevidence)*round(exp(logSamples))
  integerPart = trunc(multiplicity)
  decimalPart = multiplicity - integerPart
  
  for (mm in 1:numNestedSamples) {
    if(runif(1) <= decimalPart[mm]) integerPart[mm] = integerPart[mm] +1
  }
  
  equallyWeightedSamples = numeric()## still not efficient
  
  ## Remove logWeights from posterior as now each point should be equally weighted
  for (i in 1:numNestedSamples) {
    equallyWeightedSamples = c(equallyWeightedSamples, rep(posterior[i,-1], times = integerPart[i]))## Make a test for this as it could quite likely go wrong
  }
  
  equallyWeightedSamples = matrix(equallyWeightedSamples,
                                  ncol = ncol(posterior) - 1L, byrow=TRUE)
  
  colnames(equallyWeightedSamples) = c("logLikelihood", paste0("parameter",seq(2, ncol(equallyWeightedSamples)) -1L))
  
  return(equallyWeightedSamples)
  ### equallyWeightedSamples: A matrix of equally weighted posterior
  ### samples, with their log likelihood and parameter values in
  ### columns.
}

getEqualSamples = structure(function(
  ### Return n equally weighted posterior samples
  posterior,
  ### Matrix from the output of nested sampling. Should have log weights
  ### in column 1, log likelihoods in column 2 and then the parameter
  ### values in the remaining column(s).
  n=Inf
  ### Number of samples from the posterior required. If infinity (the deafult) this
  ### will return the maximum number of equally weighted samples
  ### generated from the posterior it can. Likewise if n is greater than
  ### this maximum number of samples.
) {
  eqSamps = .staircaseSampling(posterior)
  if (n > nrow(eqSamps) | is.infinite(n)) {
    return(eqSamps)
  } else {
    chosenRows = sample(nrow(eqSamps), size = n, replace = FALSE)
    chosenSamples = eqSamps[chosenRows, ]
    return(chosenSamples)
  }
  ### A set of equally weighted samples from the inferred posterior
  ### distribution.
}, ex=function() {

  mu <- 1
  sigma <- 1
  data <- rnorm(100, mu, sigma)
  
  transform <- function(params) {
    tParams = numeric(length=length(params))
    tParams[1] = GaussianPrior(params[1], mu, sigma)
    tParams[2] = UniformPrior(params[2], 0, 2 * sigma)
    return(tParams)
  }

  llf <- function(params) {
    tParams = transform(params)
    mean = tParams[1]
    sigma = tParams[2]
    n <- length(data)
    ll <- -(n/2) * log(2*pi) - (n/2) * log(sigma**2) - (1/(2*sigma**2)) * sum((data-mean)**2)
    return(ll)
  }
  
  prior.size <- 25
  tol <- 0.5

  ns.results <- nestedSampling(llf, 2, prior.size, transform, tolerance=tol)

  getEqualSamples(ns.results$posterior)
})


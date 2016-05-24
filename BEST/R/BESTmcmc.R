
### Parallel processing version using jagsUI::jags.basic

BESTmcmc <-
function( y1, y2=NULL, priors=NULL, doPriorsOnly=FALSE,
    numSavedSteps=1e5, thinSteps=1, burnInSteps = 1000,
    verbose=TRUE, rnd.seed=NULL, parallel=NULL) {
  # This function generates an MCMC sample from the posterior distribution.
  # y1, y2 the data vectors; y2=NULL if only one group.
  # priors is a list specifying priors to use.
  # verbose=FALSE suppresses output to the R Console.
  # rnd.seed is passed to each of the chains, with a different pseudorandom
  #   number generator for each.
  # Returns a data frame, not a matrix, with class 'BEST',
  #   with attributes Rhat, n.eff, a list with the original data, and the priors.
  #------------------------------------------------------------------------------

  if(doPriorsOnly && verbose)
    cat("Warning: The output shows the prior distributions,
      NOT the posterior distributions for your data.\n")
  # Parallel processing check
  nCores <- detectCores()
  if(!is.null(parallel) && parallel && nCores < 4)  {
    if(verbose)
      warning("Not enough cores for parallel processing, running chains sequentially.")
    parallel <- FALSE
  }
  if(is.null(parallel))
    parallel <- nCores > 3

  # Data checks
  y <- c( y1 , y2 ) # combine data into one vector
  if(!all(is.finite(y)))
    stop("The input data include NA or Inf.")
  if(length(unique(y)) < 2 &&      # sd(y) will be 0 or NA; ok if priors specified.
        (is.null(priors) ||
          is.null(priors$muSD) || 
          is.null(priors$sigmaMode) || 
          is.null(priors$sigmaSD)))
  stop("If priors are not specified, data must include at least 2 (non-equal) values.")
  
  # Prior checks:
  if(!is.null(priors))  {
    if(!is.list(priors)) {
      if(is.numeric(priors)) {
        stop("'priors' is now the 3rd argument; it must be a list (or NULL).")
      } else {
        stop("'priors' must be a list (or NULL).")
      }
    }
    nameOK <- names(priors) %in%
          c("muM", "muSD", "sigmaMode", "sigmaSD", "nuMean", "nuSD")
    if(!all(nameOK))
      stop("Invalid items in prior specification: ",
          paste(sQuote(names(priors)[!nameOK]), collapse=", "))
    if(!all(sapply(priors, is.numeric)))
      stop("All items in 'priors' must be numeric.")
    if(!is.null(priors$muSD) && priors$muSD <= 0)
      stop("muSD must be > 0")
  }
  if(is.null(rnd.seed))
    rnd.seed <- floor(runif(1,1,10000))

  # THE PRIORS
  if(is.null(priors)) {   # use the old prior specification
    dataForJAGS <- list(
      muM = mean(y) ,
      muP = 0.000001 * 1/sd(y)^2 ,
      sigmaLow = sd(y) / 1000 ,
      sigmaHigh = sd(y) * 1000
    )
  } else {    # use gamma priors
    priors0 <- list(  # default priors
      muM = mean(y) ,
      muSD = sd(y)*5 ,
      sigmaMode = sd(y),
      sigmaSD = sd(y)*5,
      nuMean = 30,
      nuSD = 30 )
    priors0 <- modifyList(priors0, priors)  # user's priors take prior-ity (duh!!)
    # Convert to Shape/Rate
    sigmaShRa <- gammaShRaFromModeSD(mode=priors0$sigmaMode, sd=priors0$sigmaSD)
    nuShRa <- gammaShRaFromMeanSD(mean=priors0$nuMean, sd=priors0$nuSD)
    dataForJAGS <- list(
      muM = priors0$muM,
      muP = 1/priors0$muSD^2,  # convert SD to precision
      Sh = sigmaShRa$shape,
      Ra = sigmaShRa$rate)
    if(!is.null(y2))  {   # all the above must be vectors of length 2
      fixPrior <- function(x) {
        if(length(x) < 2)
          x <- rep(x, 2)
        return(x)
      }
      dataForJAGS <- lapply(dataForJAGS, fixPrior)
    }
    dataForJAGS$ShNu <- nuShRa$shape
    dataForJAGS$RaNu <- nuShRa$rate
  }

  # THE MODEL.
  modelFile <- file.path(tempdir(), "BESTmodel.txt")
  if(is.null(priors)) {  # use old broad priors
    if(is.null(y2)) {
      modelString = "
      model {
        for ( i in 1:Ntotal ) {
          y[i] ~ dt( mu , tau , nu )
        }
        mu ~ dnorm( muM , muP )
        tau <- 1/pow( sigma , 2 )
        sigma ~ dunif( sigmaLow , sigmaHigh )
        nu <- nuMinusOne+1
        nuMinusOne ~ dexp(1/29)
      }
      " # close quote for modelString
    } else {
      modelString <- "
      model {
        for ( i in 1:Ntotal ) {
          y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )
        }
        for ( j in 1:2 ) {
          mu[j] ~ dnorm( muM , muP )
          tau[j] <- 1/pow( sigma[j] , 2 )
          sigma[j] ~ dunif( sigmaLow , sigmaHigh )
        }
        nu <- nuMinusOne+1
        nuMinusOne ~ dexp(1/29)
      }
      " # close quote for modelString
    }
  } else {    # use gamma priors
    if(is.null(y2)) {
      modelString = "
      model {
        for ( i in 1:Ntotal ) {
          y[i] ~ dt( mu , tau , nu )
        }
        mu ~ dnorm( muM[1] , muP[1] )
        tau <- 1/pow( sigma , 2 )
         sigma ~ dgamma( Sh[1] , Ra[1] )
        nu ~ dgamma( ShNu , RaNu ) # prior for nu
      }
      " # close quote for modelString
    } else {
      modelString <- "
      model {
        for ( i in 1:Ntotal ) {
          y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )
        }
        for ( j in 1:2 ) {
          mu[j] ~ dnorm( muM[j] , muP[j] )
          tau[j] <- 1/pow( sigma[j] , 2 )
          sigma[j] ~ dgamma( Sh[j] , Ra[j] )
        }
        nu ~ dgamma( ShNu , RaNu ) # prior for nu
      }
      " # close quote for modelString
    }
  }
  # Write out modelString to a text file
  writeLines( modelString , con=modelFile )

  #------------------------------------------------------------------------------
  # THE DATA.
  # dataForJAGS already has the priors, add the data:
  if(!doPriorsOnly)
    dataForJAGS$y <- y
  dataForJAGS$Ntotal <- length(y)
  if(!is.null(y2)) # create group membership code
    dataForJAGS$x <- c( rep(1,length(y1)) , rep(2,length(y2)) )

  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  if(is.null(y2)) {
    mu = mean(y1)
    sigma = sd(y1)
  } else {
    mu = c( mean(y1) , mean(y2) )
    sigma = c( sd(y1) , sd(y2) )
  }
  # Regarding initial values in next line: (1) sigma will tend to be too big if
  # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  # initial values keep the burn-in period moderate.
  
  initsList0 <- list(mu=mu, sigma=sigma, .RNG.seed=rnd.seed)
  if(is.null(priors)) {
    initsList0$nuMinusOne <- 4
  } else {
    initsList0$nu <- 5
  }
  initsList <- list(
                c(initsList0, .RNG.name="base::Wichmann-Hill"),
                c(initsList0, .RNG.name="base::Marsaglia-Multicarry"),
                c(initsList0, .RNG.name="base::Super-Duper") )
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  codaSamples <- jags.basic(
    data = dataForJAGS,
    inits = initsList,
    parameters.to.save = c( "mu" , "sigma" , "nu" ),     # The parameters to be monitored
    model.file = modelFile,
    n.chains = 3,    # Do not change this without also changing inits.
    n.adapt = 500,
    n.iter = ceiling( ( numSavedSteps * thinSteps) / 3  + burnInSteps ),
    n.burnin = burnInSteps,
    n.thin = thinSteps,
    modules = NULL,
    parallel = parallel,
    DIC = FALSE,
    seed = rnd.seed,
    verbose = verbose)
  #------------------------------------------------------------------------------
  mcmcChain = as.matrix( codaSamples )
  if(dim(mcmcChain)[2] == 5 &&
        all(colnames(mcmcChain) == c("mu[1]", "mu[2]", "nu", "sigma[1]", "sigma[2]")))
    colnames(mcmcChain) <- c("mu1", "mu2", "nu", "sigma1", "sigma2")
  mcmcDF <- as.data.frame(mcmcChain)
  class(mcmcDF) <- c("BEST", class(mcmcDF))
  attr(mcmcDF, "call") <- match.call()
  attr(mcmcDF, "Rhat") <- gelman.diag(codaSamples)$psrf[, 1]
  attr(mcmcDF, "n.eff") <- effectiveSize(codaSamples)
  attr(mcmcDF, "data") <- list(y1 = y1, y2 = y2)
  attr(mcmcDF, "doPriorsOnly") <- doPriorsOnly
  if(!is.null(priors))
    attr(mcmcDF, "priors") <- priors0

  return( mcmcDF )
}

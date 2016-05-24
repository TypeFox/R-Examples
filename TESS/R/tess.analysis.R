################################################################################
#
# tess.analysis.R
#
# Copyright (c) 2012- Sebastian Hoehna
#
# This file is part of TESS.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# TESS is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  TESS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with TESS; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA
#
################################################################################

################################################################################
#
# @brief Running an analysis on a given tree and estimating the diversification
#        rates including rate-shifts and mass-extinction events.
#
# @date Last modified: 2014-10-05
# @author Sebastian Hoehna
# @version 2.0
# @since 2012-09-22, version 2.0
#
# @param    tree                                          phylo         the tree
# @param    initialSpeciationRate                         vector        The initial value of the speciation rate when the MCMC is started.
# @param    initialExtinctionRate                         vector        The initial value of the extinction rate when the MCMC is started.
# @param    speciationRatePriorMean                       scalar        mean parameter for the lognormal prior on lambda
# @param    speciationRatePriorStDev                      scalar        variance parameter for the lognormal prior on lambda
# @param    extinctionRatePriorMean                       scalar        mean parameter for the lognormal prior on mu
# @param    extinctionRatePriorStDev                      scalar        variance parameter for the lognormal prior on mu
# @param    initialSpeciationRateChangeTime               vector        The initial value of the time points when speciation rate-shifts occur.
# @param    initialExtinctionRateChangeTime               vector        The initial value of the time points when speciation rate-shifts occur.
# @param    estimateNumberRateChanges                     boolean       Do we estimate the number of rate shifts?
# @param    numExpectedRateChanges                        scalar        Expected number of rate changes which follow a Poisson process.
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    missingSpecies                                vector        The number of species missed which originated in a given time interval (empirical taxon sampling)
# @param    timesMissingSpecies                           vector        The times intervals of the missing species (empirical taxon sampling)
# @param    tInitialMassExtinction                        vector        The initial value of the vector of times of the mass-extinction events.
# @param    pInitialMassExtinction                        vector        The initial value of the vector of survival probabilities of the mass-extinction events.
# @param    pMassExtinctionPriorShape1                    scalar        The alpha (first shape) parameter of the Beta prior distribution for the survival probability of a mass-extinction event.
# @param    pMassExtinctionPriorShape2                    scalar        The beta (second shape) parameter of the Beta prior distribution for the survival probability of a mass-extinction event.
# @param    estimateMassExtinctionTimes                   boolean       Do we estimate the times of mass-extinction events?
# @param    numExpectedMassExtinctions                    scalar        Expected number of mass-extinction events which follow a Poisson process.
# @param    estimateNumberMassExtinctions                 boolean       Do we estimate the number of mass-extinction events?
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @param    MAX_ITERATIONS                                scalar        The maximum number of iteration of the MCMC.
# @param    THINNING                                      scalar        The frequency how often samples are recorded during the MCMC.
# @param    OPTIMIZATION_FREQUENCY                        scalar        The frequency how often the MCMC moves are optimized
# @param    CONVERGENCE_FREQUENCY                         scalar        The frequency how often we check for convergence?
# @param    MAX_TIME                                      scalar        The maximum time the MCMC is allowed to run in seconds.
# @param    MIN_ESS                                       scalar        The minimum number of effective samples (ESS) to assume convergence.
# @param    ADAPTIVE                                      scalar        Do we use auto-tuning of the MCMC moves?
# @param    dir                                           string        The subdirectory in which the output will be stored.
# @param    priorOnly                                     boolean       Do we sample from the prior only?
# @param    verbose                                       boolean       Do we want to print the progress of the MCMC?
#
################################################################################

tess.analysis <- function( tree,
                           initialSpeciationRate,
                           initialExtinctionRate,
                           empiricalHyperPriors = TRUE,
                           empiricalHyperPriorInflation = 10.0,
                           empiricalHyperPriorForm = c("lognormal","normal","gamma"),
                           speciationRatePriorMean = 0.0,
                           speciationRatePriorStDev = 1.0,
                           extinctionRatePriorMean = 0.0,
                           extinctionRatePriorStDev = 1.0,
                           initialSpeciationRateChangeTime = c(),
                           initialExtinctionRateChangeTime = c(),
                           estimateNumberRateChanges = TRUE,
                           numExpectedRateChanges = 2,
                           samplingProbability = 1,
                           missingSpecies = c(),
                           timesMissingSpecies = c(),
                           tInitialMassExtinction = c(),
                           pInitialMassExtinction = c(),
                           pMassExtinctionPriorShape1 = 5,
                           pMassExtinctionPriorShape2 = 95,
                           estimateMassExtinctionTimes = TRUE,
                           numExpectedMassExtinctions = 2,
                           estimateNumberMassExtinctions = TRUE,
                           MRCA = TRUE,
                           CONDITION = "survival",
                           BURNIN = 10000,
                           MAX_ITERATIONS = 200000,
                           THINNING = 100,
                           OPTIMIZATION_FREQUENCY = 500,
                           CONVERGENCE_FREQUENCY = 1000,
                           MAX_TIME = Inf, MIN_ESS = 500,
                           ADAPTIVE = TRUE,
                           dir = "" ,
                           priorOnly = FALSE,
                           verbose = TRUE) {

  orgDir <- getwd()
  if ( dir != "" ) {
    if (file.exists(dir)){
      setwd(file.path(dir))
    } else {
      dir.create(file.path(dir))
      setwd(file.path(dir))
    }
  }

  if ( BURNIN > 0 & BURNIN < 1 ) {
    burn <- MAX_ITERATIONS * BURNIN
  } else {
    burn <- BURNIN
  }

  MAX_TRIES <- 1000

  times <- branching.times(tree)

  if ( !is.ultrametric(tree) ) {
    stop("The likelihood function is only defined for ultrametric trees!")
  }

  write.nexus(tree,file="input_tree.tre")

  if ( empiricalHyperPriors == TRUE ) {

    if( verbose ) {
      cat("Estimating empirical hyper-parameters.\n\n")
    }

    empirical.prior.likelihood <- function(params) {
      # We use the parameters as diversification rate and turnover rate.
      # Thus we need to transform first
#       b <- params[1] + params[2]
#       d <- params[2]
      b <- params[1]/(1-params[2])
      d <- params[2]*b

      lnl <- tess.likelihood(times, b, d, missingSpecies = missingSpecies, timesMissingSpecies = timesMissingSpecies, samplingStrategy = "uniform", samplingProbability= samplingProbability, MRCA = MRCA,CONDITION = CONDITION,log=TRUE)
      return(lnl)
    }

    prior.diversification <- function(x) { dexp(x,rate=max(times)/log(length(times)),log=TRUE) }
    prior.turnover <- function(x) { dbeta(x,shape1=1.5,shape2=3,log=TRUE) }
    priors <- c("diversification"=prior.diversification,"turnover"=prior.turnover)

    tries <- 1
    repeat {
        starting.values <- runif(2,0,1)
        starting.lnl <- empirical.prior.likelihood(starting.values)
        if ( is.finite(starting.lnl) || tries >= MAX_TRIES ) {
            break
        }
        tries <- tries + 1
    }

    if ( tries >= MAX_TRIES ) {
        stop("Could not find valid starting values after ",tries," attemps.\n")
    }

    samples <- tess.mcmc( empirical.prior.likelihood,priors,starting.values,logTransforms=c(TRUE,TRUE),delta=c(0.1,0.1),iterations=20000,burnin=2000,verbose=verbose)

    samples.lambda <- samples[,1]/(1-samples[,2])
    m.lambda <- mean(samples.lambda)
    v.lambda <- empiricalHyperPriorInflation * var(samples.lambda)

    samples.mu <- samples.lambda*samples[,2]
    m.mu <- mean(samples.mu)
    v.mu <- empiricalHyperPriorInflation * var(samples.mu)

    lambda.hyperprior.fits <- c("lognormal"=Inf,"normal"=Inf,"gamma"=Inf)
    mu.hyperprior.fits <- c("lognormal"=Inf,"normal"=Inf,"gamma"=Inf)

    if ( "lognormal" %in% empiricalHyperPriorForm ) {
      lambda.hyperprior.fits["lognormal"] <- suppressWarnings(ks.test(samples.lambda,'plnorm',meanlog=log((m.lambda^2)/sqrt(v.lambda+m.lambda^2)),sdlog=sqrt( log(1+v.lambda/(m.lambda^2)))))$statistic
      mu.hyperprior.fits["lognormal"] <- suppressWarnings(ks.test(samples.mu,'plnorm',meanlog=log((m.mu^2)/sqrt(v.mu+m.mu^2)),sdlog=sqrt( log(1+v.mu/(m.mu^2)))))$statistic
    }

    if ( "normal" %in% empiricalHyperPriorForm ) {
      lambda.hyperprior.fits["normal"] <- suppressWarnings(ks.test(samples.lambda,'pnorm',mean=m.lambda,sd=sqrt(v.lambda)))$statistic
      mu.hyperprior.fits["normal"] <- suppressWarnings(ks.test(samples.mu,'pnorm',mean=m.mu,sd=sqrt(v.mu)))$statistic
    }

    if ( "gamma" %in% empiricalHyperPriorForm ) {
      lambda.hyperprior.fits["gamma"] <- suppressWarnings(ks.test(samples.lambda,'pgamma',shape=(m.lambda^2)/v.lambda,rate=m.lambda/v.lambda))$statistic
      mu.hyperprior.fits["gamma"] <- suppressWarnings(ks.test(samples.mu,'pgamma',shape=(m.mu^2)/v.mu,rate=m.mu/v.mu))$statistic
    }

    lambda.hyper.prior.form <- names(which.min(lambda.hyperprior.fits[empiricalHyperPriorForm]))
    mu.hyper.prior.form <- names(which.min(mu.hyperprior.fits[empiricalHyperPriorForm]))

    if ( lambda.hyper.prior.form == "lognormal" ) {
      speciationRatePriorMean <- log((m.lambda^2)/sqrt(v.lambda+m.lambda^2))
      speciationRatePriorStDev <- sqrt( log(1+v.lambda/(m.lambda^2)) )
    } else if ( lambda.hyper.prior.form == "normal" ) {
      speciationRatePriorMean <- m.lambda
      speciationRatePriorStDev <- sqrt(v.lambda)
    } else if ( lambda.hyper.prior.form == "gamma" ) {
      speciationRatePriorMean <- (m.lambda^2)/v.lambda
      speciationRatePriorStDev <- m.lambda/v.lambda
    }

    if ( mu.hyper.prior.form == "lognormal" ) {
      extinctionRatePriorMean <- log((m.mu^2)/sqrt(v.mu+m.mu^2))
      extinctionRatePriorStDev <- sqrt( log(1+v.mu/(m.mu^2)))
    } else if ( mu.hyper.prior.form == "normal" ) {
      extinctionRatePriorMean <- m.mu
      extinctionRatePriorStDev <- sqrt(v.mu)
    } else if ( mu.hyper.prior.form == "gamma" ) {
      extinctionRatePriorMean <- (m.mu^2)/v.mu
      extinctionRatePriorStDev <- m.mu/v.mu
    }

    initialSpeciationRate <- m.lambda
    initialExtinctionRate <- m.mu

    empirical.hyperprior.samples <- data.frame(1:nrow(samples),samples.lambda,samples.mu)
    write.table(empirical.hyperprior.samples[-1,],file='samplesHyperprior.txt',row.names=FALSE,col.names=c('Iteration','speciation','extinction'),sep='\t')

    summary.file <- "empiricalHyperpriorSummary.txt"
    cat("Summary of the empirical hyperprior analysis",file=summary.file)
    cat("\n============================================\n\n",file=summary.file,append=TRUE)
    cat("Specation rate hyperprior\n",file=summary.file,append=TRUE)
    cat("-------------------------\n",file=summary.file,append=TRUE)
    cat("Form:\t",lambda.hyper.prior.form,'\n',file=summary.file,append=TRUE)
    if ( lambda.hyper.prior.form == "lognormal" ) {
      cat("meanlog:\t",speciationRatePriorMean,"\n",file=summary.file,append=TRUE)
      cat("sdlog:\t\t",speciationRatePriorStDev,"\n",file=summary.file,append=TRUE)
    } else if ( lambda.hyper.prior.form == "normal" ) {
      cat("mean:\t",speciationRatePriorMean,"\n",file=summary.file,append=TRUE)
      cat("sd:\t\t",speciationRatePriorStDev,"\n",file=summary.file,append=TRUE)
    } else if ( lambda.hyper.prior.form == "gamma" ) {
      cat("shape:\t",speciationRatePriorMean,"\n",file=summary.file,append=TRUE)
      cat("rate:\t\t",speciationRatePriorStDev,"\n",file=summary.file,append=TRUE)
    }
    cat("\nExtinction rate hyperprior\n",file=summary.file,append=TRUE)
    cat("-------------------------\n",file=summary.file,append=TRUE)
    cat("Form:\t",mu.hyper.prior.form,'\n',file=summary.file,append=TRUE)
    if ( mu.hyper.prior.form == "lognormal" ) {
      cat("meanlog:\t",extinctionRatePriorMean,"\n",file=summary.file,append=TRUE)
      cat("sdlog:\t\t",extinctionRatePriorStDev,"\n",file=summary.file,append=TRUE)
    } else if ( mu.hyper.prior.form == "normal" ) {
      cat("mean:\t",extinctionRatePriorMean,"\n",file=summary.file,append=TRUE)
      cat("sd:\t\t",extinctionRatePriorStDev,"\n",file=summary.file,append=TRUE)
    } else if ( mu.hyper.prior.form == "gamma" ) {
      cat("shape:\t",extinctionRatePriorMean,"\n",file=summary.file,append=TRUE)
      cat("rate:\t\t",extinctionRatePriorStDev,"\n",file=summary.file,append=TRUE)
    }

  } else { # no empirical hyper priors

    lambda.hyper.prior.form <- mu.hyper.prior.form <- empiricalHyperPriorForm[1]

  }

  likelihood <- function(lambda,lambdaTimes,mu,muTimes,tMassExtinction,pSurvival) {

    if ( priorOnly == FALSE ) {
      if( any( lambda < 0 ) | any( mu < 0 ) ) {
        lnl <- -Inf
      } else {
        lnl <- tess.likelihood.rateshift(times, lambda, mu, rateChangeTimesLambda = lambdaTimes, rateChangeTimesMu = muTimes, massExtinctionTimes= tMassExtinction, massExtinctionSurvivalProbabilities = pSurvival, missingSpecies = missingSpecies, timesMissingSpecies = timesMissingSpecies, samplingStrategy = "uniform", samplingProbability= samplingProbability, MRCA = MRCA,CONDITION = CONDITION,log = TRUE)
      }
    } else {
      lnl <- 0
    }
    return (lnl)

  }

  prior <- function(timesLambda,valuesLambda,timesMu,valuesMu,timesMassExtinction,valuesMassExtinction) {

    if( any( valuesLambda < 0 ) | any( valuesMu < 0 ) ) {
      return (-Inf)
    }

    #####################
    ## Speciation Rate ##
    #####################

    # prior on number of changes for lambda
    k <- length(timesLambda)
    lnp <- dpois(k,numExpectedRateChanges,log=TRUE)

    if ( k > 0 ) {
      # prior on change times
      lnp <- lnp + sum( dunif(timesLambda,0,AGE,log=TRUE) )
    }

    # prior on change values
    if ( lambda.hyper.prior.form == "lognormal" ) {
      lnp <- lnp + sum( dlnorm(valuesLambda,speciationRatePriorMean,speciationRatePriorStDev,log=TRUE) )
    } else if ( lambda.hyper.prior.form == "normal" ) {
      lnp <- lnp + sum( dnorm(valuesLambda,speciationRatePriorMean,speciationRatePriorStDev,log=TRUE) )
    } else if ( lambda.hyper.prior.form == "gamma" ) {
      lnp <- lnp + sum( dgamma(valuesLambda,speciationRatePriorMean,speciationRatePriorStDev,log=TRUE) )
    }

    #####################
    ## Extinction Rate ##
    #####################

    # prior on number of changes for mu
    k <- length(timesMu)
    lnp <- lnp + dpois(k,numExpectedRateChanges,log=TRUE)

    if ( k > 0 ) {
      # prior on change times
      lnp <- lnp + sum( dunif(timesMu,0,AGE,log=TRUE) )
    }

    # prior on change values
    if ( mu.hyper.prior.form == "lognormal" ) {
      lnp <- lnp + sum( dlnorm(valuesMu,extinctionRatePriorMean,extinctionRatePriorStDev,log=TRUE) )
    } else if ( mu.hyper.prior.form == "normal" ) {
      lnp <- lnp + sum( dnorm(valuesMu,extinctionRatePriorMean,extinctionRatePriorStDev,log=TRUE) )
    } else if ( mu.hyper.prior.form == "gamma" ) {
      lnp <- lnp + sum( dgamma(valuesMu,extinctionRatePriorMean,extinctionRatePriorStDev,log=TRUE) )
    }

    #####################
    ## Mass-Extinction ##
    #####################

    # prior on number of changes for mass-extinction events
    k <- length(timesMassExtinction)
    lnp <- lnp + dpois(k,numExpectedMassExtinctions,log=TRUE)

    if ( k > 0 ) {
      # prior on change times
      lnp <- lnp + sum( dunif(timesMassExtinction,0,AGE,log=TRUE) )

      # prior on change values
      lnp <- lnp + sum( dbeta(valuesMassExtinction,pMassExtinctionPriorShape1,pMassExtinctionPriorShape2,log=TRUE) )
    }

    return (lnp)

  }

  AGE <- max( as.numeric( branching.times( tree ) ) )

  # these are the parameters we sample
  posterior <- c()
  kLambda <- c()
  kMu <- c()
  kMassExtinction <- c()
  speciationRateValues <- list()
  speciationRateChangeTimes <- list()
  extinctionRateValues <- list()
  extinctionRateChangeTimes <- list()
  survivalProbability <- list()
  massExtinctionTime <- list()

  cat("SpeciationRates\n",file="SpeciationRates.txt")
  cat("SpeciationRateChanges\n",file="SpeciationRateChanges.txt")
  cat("ExtinctionRates\n",file="ExtinctionRates.txt")
  cat("ExtinctionRateChangeTimes\n",file="ExtinctionRateChanges.txt")
  cat("SurvivalProbabilities\n",file="SurvivalProbabilities.txt")
  cat("MassExtinctionTimes\n",file="MassExtinctionTimes.txt")
  cat(paste("Iteration","posterior","NumSpeciation","numExtinction","numMassExtinctions",sep="\t"),"\n",file="samples_numCategories.txt")

  # the "random" starting values
  lambda <- initialSpeciationRate
  mu <- initialExtinctionRate
  lambdaChangeTimes <- initialSpeciationRateChangeTime
  muChangeTimes <- initialExtinctionRateChangeTime
  pMassExtinction <- pInitialMassExtinction
  tMassExtinction <- tInitialMassExtinction
  initialPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)
  currentPP <- initialPP

  ## check if the prior and posterior gives valid probabilities
  MAX_TRIES <- 1000
  for ( i in 1:MAX_TRIES ) {
    tmpPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)
    if ( is.finite(tmpPP) ) break
    lambdaChangeTimes <- c()
    lambda <- rlnorm(1,speciationRatePriorMean,speciationRatePriorStDev)
    muChangeTimes <- c()
    mu <- rlnorm(1,extinctionRatePriorMean,extinctionRatePriorStDev)
    tMassExtinction <- c()
    pMassExtinction <- c()
  }

  ## MCMC proposal sizes
  deltaLambdaTime <- 0.1
  lambdaTimeAccepted <- 0
  lambdaTimeTried <- 0
  deltaLambdaValue <- 0.1
  lambdaValueAccepted <- 0
  lambdaValueTried <- 0
  deltaMuTime <- 0.1
  muTimeAccepted <- 0
  muTimeTried <- 0
  deltaMuValue <- 0.1
  muValueAccepted <- 0
  muValueTried <- 0
  deltaMassExtinctionTime <- 0.1
  massExtinctionTimeAccepted <- 0
  massExtinctionTimeTried <- 0
  deltaMassExtinctionValue <- 0.1
  massExtinctionValueAccepted <- 0
  massExtinctionValueTried <- 0

  finished <- FALSE
  SAMPLE <- FALSE
  min.ess <- 0
  startTime <- Sys.time()
  i <- 1

  if ( verbose ) {
    cat("\nPerforming CoMET analysis.\n\n")
    if( burn > 0 ) {
      cat("Burning-in the chain ...\n")
    } else {
      cat("Running the chain ... \n")
    }
    cat("0--------25--------50--------75--------100\n")
    bar <- txtProgressBar(style=1,width=42)
  }

  while ( !finished ) {

    # change the time of an event
    if ( length(lambdaChangeTimes) > 0 ) {

      # compute index of value that will be changed
      idx <- floor( runif(1,0,length(lambdaChangeTimes)) ) + 1

      # store current value
      t <- lambdaChangeTimes[idx]

      # compute current posterior
      oldPP <- currentPP

      # propose new value
      tPrime <- rnorm(1,t,deltaLambdaTime)
      lambdaChangeTimes[idx] <- tPrime

      # compute new posterior
      if ( tPrime > 0 && tPrime < AGE ) {
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)
      } else {
        newPP <- -Inf
      }

      # accept/reject
      if ( is.finite(newPP - oldPP) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP ) ) {
        # reject
        lambdaChangeTimes[idx] <- t
      } else {
        currentPP <- newPP
        lambdaTimeAccepted <- lambdaTimeAccepted + 1
      }
      lambdaTimeTried <- lambdaTimeTried + 1
    }

    # change the value of an event

    # compute index of value that will be changed
    idx <- floor( runif(1,0,length(lambda)) ) + 1

    # store current value
    v <- lambda[idx]

    # compute current posterior
    oldPP <- currentPP

    # propose new value
    u <- rnorm(1,0,deltaLambdaValue)
    scalingFactor <- exp( u )
    vPrime <- v * scalingFactor
    lambda[idx] <- vPrime

    # compute the Hastings ratio
    lnHastingsratio <- log( scalingFactor )

    # compute new posterior
    newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

    # accept/reject
    if ( is.finite(newPP - oldPP + lnHastingsratio) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + lnHastingsratio ) ) {
      # reject
      lambda[idx] <- v
    } else {
      currentPP <- newPP
      lambdaValueAccepted <- lambdaValueAccepted + 1
    }
    lambdaValueTried <- lambdaValueTried + 1

    if ( estimateNumberRateChanges == TRUE ) {
      # We need to randomly pick a birth or death move
      # Otherwise we might give birth and die every time
      u <- runif(1,0,1)
      if ( u > 0.5 ) {
        # birth move

        # compute current posterior
        oldPP <- currentPP

        # randomly pick a new time
        t <- runif(1,0,AGE)

        # randomly pick a new value
        if ( lambda.hyper.prior.form == "lognormal" ) {
          v <- rlnorm(1,speciationRatePriorMean,speciationRatePriorStDev)
        } else if ( lambda.hyper.prior.form == "normal" ) {
          v <- rnorm(1,speciationRatePriorMean,speciationRatePriorStDev)
        } else if ( lambda.hyper.prior.form == "gamma" ) {
          v <- rgamma(1,speciationRatePriorMean,speciationRatePriorStDev)
        }

        # construct the new parameters
        lambdaChangeTimes[length(lambdaChangeTimes)+1] <- t
        lambda[length(lambda)+1] <- v

        # compute new posterior
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

        # compute proposal ratio
        kPrime <- length(lambdaChangeTimes)
        if ( lambda.hyper.prior.form == "lognormal" ) {
          pr <- 1 / ( dunif(t,0,AGE) * dlnorm(v,speciationRatePriorMean,speciationRatePriorStDev) )
        } else if ( lambda.hyper.prior.form == "normal" ) {
          pr <- 1 / ( dunif(t,0,AGE) * dnorm(v,speciationRatePriorMean,speciationRatePriorStDev) )
        } else if ( lambda.hyper.prior.form == "gamma" ) {
          pr <- 1 / ( dunif(t,0,AGE) * dgamma(v,speciationRatePriorMean,speciationRatePriorStDev) )
        }

        # accept/reject
        if ( is.finite(newPP - oldPP + log(pr)) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + log( pr ) ) ) {
          # reject
          lambdaChangeTimes <- lambdaChangeTimes[-length(lambdaChangeTimes)]
          lambda <- lambda[-length(lambda)]
        } else {
          currentPP <- newPP
        }
      } else {

        # death move
        if ( length(lambdaChangeTimes) > 0 ) {

          # compute current posterior
          oldPP <- currentPP

          # randomly pick an index
          idx <- floor( runif(1,0,length(lambdaChangeTimes)) ) + 1

          # store the current value
          t <- lambdaChangeTimes
          v <- lambda

          # construct the new parameters
          lambdaChangeTimes <- lambdaChangeTimes[-idx]
          lambda <- lambda[-(idx+1)]

          # compute new posterior
          newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

          # compute proposal ratio
          kPrime <- length(lambdaChangeTimes) + 1
          if ( lambda.hyper.prior.form == "lognormal" ) {
            pr2 <-  dunif(t[idx],0,AGE) * dlnorm(v[idx+1],speciationRatePriorMean,speciationRatePriorStDev)
          } else if ( lambda.hyper.prior.form == "normal" ) {
            pr2 <-  dunif(t[idx],0,AGE) * dnorm(v[idx+1],speciationRatePriorMean,speciationRatePriorStDev)
          } else if ( lambda.hyper.prior.form == "gamma" ) {
            pr2 <-  dunif(t[idx],0,AGE) * dgamma(v[idx+1],speciationRatePriorMean,speciationRatePriorStDev)
          }

          # accept/reject
          if ( is.finite(newPP - oldPP + log(pr2)) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + log( pr2 ) ) ) {
            # reject
            lambdaChangeTimes <- t
            lambda <- v
          } else {
            currentPP <- newPP
          }
        }
      }
    }

    #####################
    ## Extinction Rate ##
    #####################

    # change the time of an event
    if ( length(muChangeTimes) > 0 ) {

      # compute index of value that will be changed
      idx <- floor( runif(1,0,length(muChangeTimes)) ) + 1

      # store current value
      t <- muChangeTimes[idx]

      # compute current posterior
      oldPP <- currentPP

      # propose new value
      tPrime <- rnorm(1, t, deltaMuTime)
      muChangeTimes[idx] <- tPrime

      # compute new posterior
      if ( tPrime > 0 && tPrime < AGE ) {
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)
      } else {
        newPP <- -Inf
      }

      # accept/reject
      if ( is.finite(newPP - oldPP) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP ) ) {
        # reject
        muChangeTimes[idx] <- t
      } else {
        currentPP <- newPP
        muTimeAccepted <- muTimeAccepted + 1
      }
      muTimeTried <- muTimeTried + 1
    }

    # change the value of an event

    # compute index of value that will be changed
    idx <- floor( runif(1,0,length(mu)) ) + 1

    # store current value
    v <- mu[idx]

    # compute current posterior
    oldPP <- currentPP

    # propose new value
    u <- rnorm(1,0,deltaMuValue)
    scalingFactor <- exp( u )
    vPrime <- v * scalingFactor
    mu[idx] <- vPrime

    # compute the Hastings ratio
    lnHastingsratio <- log( scalingFactor )

    # compute new posterior
    newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

    # accept/reject
    if ( is.finite(newPP - oldPP + lnHastingsratio) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + lnHastingsratio ) ) {
      # reject
      mu[idx] <- v
    } else {
      currentPP <- newPP
      muValueAccepted <- muValueAccepted + 1
    }
    muValueTried <- muValueTried + 1

    if ( estimateNumberRateChanges == TRUE ) {
      # We need to randomly pick a birth or death move
      # Otherwise we might give birth and die every time

      u <- runif(1,0,1)
      if ( u > 0.5 ) {
        # birth move

        # compute current posterior
        oldPP <- currentPP

        # randomly pick a new time
        t <- runif(1,0,AGE)

        # randomly pick a new value
        if ( mu.hyper.prior.form == "lognormal" ) {
          v <- rlnorm(1,extinctionRatePriorMean,extinctionRatePriorStDev)
        } else if ( mu.hyper.prior.form == "normal" ) {
          v <- rnorm(1,extinctionRatePriorMean,extinctionRatePriorStDev)
        } else if ( mu.hyper.prior.form == "gamma" ) {
          v <- rgamma(1,extinctionRatePriorMean,extinctionRatePriorStDev)
        }

        # construct the new parameters
        muChangeTimes[length(muChangeTimes)+1] <- t
        mu[length(mu)+1] <- v

        # compute new posterior
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

        # compute proposal ratio
        kPrime <- length(muChangeTimes)
        if ( mu.hyper.prior.form == "lognormal" ) {
          pr <- 1 / ( dunif(t,0,AGE) * dlnorm(v,extinctionRatePriorMean,extinctionRatePriorStDev) )
        } else if ( mu.hyper.prior.form == "normal" ) {
          pr <- 1 / ( dunif(t,0,AGE) * dnorm(v,extinctionRatePriorMean,extinctionRatePriorStDev) )
        } else if ( mu.hyper.prior.form == "gamma" ) {
          pr <- 1 / ( dunif(t,0,AGE) * dgamma(v,extinctionRatePriorMean,extinctionRatePriorStDev) )
        }

        # accept/reject
        if ( is.finite(newPP - oldPP + log(pr)) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + log( pr ) ) ) {
          # reject
          muChangeTimes <- muChangeTimes[-length(muChangeTimes)]
          mu <- mu[-length(mu)]
        } else {
          currentPP <- newPP
        }
      } else {

        # death move
        if ( length(muChangeTimes) > 0 ) {

          # compute current posterior
          oldPP <- currentPP

          # randomly pick an index
          idx <- floor( runif(1,0,length(muChangeTimes)) ) + 1

          # store the current value
          t <- muChangeTimes
          v <- mu

          # construct the new parameters
          muChangeTimes <- muChangeTimes[-idx]
          mu <- mu[-(idx+1)]

          # compute new posterior
          newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

          # compute proposal ratio
          kPrime <- length(muChangeTimes) + 1
          if ( mu.hyper.prior.form == "lognormal" ) {
            pr2 <-  dunif(t[idx],0,AGE) * dlnorm(v[idx+1],extinctionRatePriorMean,extinctionRatePriorStDev)
          } else if ( mu.hyper.prior.form == "normal" ) {
            pr2 <-  dunif(t[idx],0,AGE) * dnorm(v[idx+1],extinctionRatePriorMean,extinctionRatePriorStDev)
          } else if ( mu.hyper.prior.form == "gamma" ) {
            pr2 <-  dunif(t[idx],0,AGE) * dgamma(v[idx+1],extinctionRatePriorMean,extinctionRatePriorStDev)
          }

          # accept/reject
          if ( is.finite(newPP - oldPP + log(pr2)) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + log( pr2 ) ) ) {
            # reject
            muChangeTimes <- t
            mu <- v
          } else {
            currentPP <- newPP
          }
        }
      }
    }

    #####################
    ## Mass-Extinction ##
    #####################

    # change the time of an event
    if ( estimateMassExtinctionTimes == TRUE && length(tMassExtinction) > 0 ) {

      # compute index of value that will be changed
      idx <- floor( runif(1,0,length(tMassExtinction)) ) + 1

      # store current value
      t <- tMassExtinction[idx]

      # compute current posterior
      oldPP <- currentPP

      # propose new value
      tPrime <- rnorm(1, t, deltaMassExtinctionTime)
      tMassExtinction[idx] <- tPrime

      # compute new posterior
      if ( tPrime > 0 && tPrime < AGE ) {
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)
      } else {
        newPP <- -Inf
      }

      # accept/reject
      if ( is.finite(newPP - oldPP) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP ) ) {
        # reject
        tMassExtinction[idx] <- t
      } else {
        currentPP <- newPP
        massExtinctionTimeAccepted <- massExtinctionTimeAccepted + 1
      }
      massExtinctionTimeTried <- massExtinctionTimeTried + 1
    }

    # change the value of an event
    if ( length(pMassExtinction) > 0 ) {
      # compute index of value that will be changed
      idx <- floor( runif(1,0,length(pMassExtinction)) ) + 1

      # store current value
      v <- pMassExtinction[idx]

      # compute current posterior
      oldPP <- currentPP

      # propose new value
      vPrime <- rnorm(1, v, deltaMassExtinctionValue)
      pMassExtinction[idx] <- vPrime

      # compute new posterior
      if ( vPrime > 0 && vPrime < 1 ) {
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)
      } else {
        newPP <- -Inf
      }

      # accept/reject
      if ( is.finite(newPP - oldPP) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP ) ) {
        # reject
        pMassExtinction[idx] <- v
      } else {
        currentPP <- newPP
        massExtinctionValueAccepted <- massExtinctionValueAccepted + 1
      }
      massExtinctionValueTried <- massExtinctionValueTried + 1
    }

    if ( estimateNumberMassExtinctions == TRUE ) {
      # We need to randomly pick a birth or death move
      # Otherwise we might give birth and die every time

      u <- runif(1,0,1)
      if ( u > 0.5 ) {
        # birth move

        # compute current posterior
        oldPP <- currentPP

        # randomly pick a new time
        t <- runif(1,0,AGE)

        # randomly pick a new value
        v <- rbeta(1,pMassExtinctionPriorShape1,pMassExtinctionPriorShape2)

        # construct the new parameters
        tMassExtinction[length(tMassExtinction)+1] <- t
        pMassExtinction[length(pMassExtinction)+1] <- v

        # compute new posterior
        newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

        # compute proposal ratio
        kPrime <- length(pMassExtinction)
        pr <- 1 / ( dunif(t,0,AGE) * dbeta(v,pMassExtinctionPriorShape1,pMassExtinctionPriorShape2) )

        # accept/reject
        if ( is.finite(newPP - oldPP + log(pr)) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + log( pr ) ) ) {
          # reject
          tMassExtinction <- tMassExtinction[-length(tMassExtinction)]
          pMassExtinction <- pMassExtinction[-length(pMassExtinction)]
        } else {
          currentPP <- newPP
        }
      } else {

        # death move
        if ( length(tMassExtinction) > 0 ) {

          # compute current posterior
          oldPP <- currentPP

          # randomly pick an index
          idx <- floor( runif(1,0,length(tMassExtinction)) ) + 1

          # store the current value
          t <- tMassExtinction
          v <- pMassExtinction

          # construct the new parameters
          tMassExtinction <- tMassExtinction[-idx]
          pMassExtinction <- pMassExtinction[-idx]

          # compute new posterior
          newPP <- prior(lambdaChangeTimes,lambda,muChangeTimes,mu,tMassExtinction,pMassExtinction) + likelihood(lambda,lambdaChangeTimes,mu,muChangeTimes,tMassExtinction,pMassExtinction)

          # compute proposal ratio
          kPrime <- length(tMassExtinction) + 1
          pr2 <-  dunif(t[idx],0,AGE) * dbeta(v[idx],pMassExtinctionPriorShape1,pMassExtinctionPriorShape2)

          # accept/reject
          if ( is.finite(newPP - oldPP + log(pr2)) == FALSE || log( runif(1,0,1) ) > ( newPP - oldPP + log( pr2 ) ) ) {
            # reject
            tMassExtinction <- t
            pMassExtinction <- v
          } else {
            currentPP <- newPP
          }
        }
      }
    }

    #################
    ## Auto-Tuning ##
    #################

    # if we use adaptive MCMC we might have to modify our parameters
    if ( ADAPTIVE == TRUE & SAMPLE == FALSE ) {
      if ( i %% OPTIMIZATION_FREQUENCY == 0 ) {

        # delta for time of lambda change events
        if ( lambdaTimeTried > 0 ) {
          rate <- lambdaTimeAccepted / lambdaTimeTried
          if ( rate > 0.44 ) {
            deltaLambdaTime <- deltaLambdaTime * (1.0 + ((rate-0.44)/0.56) )
          } else {
            deltaLambdaTime <- deltaLambdaTime / (2.0 - rate/0.44 )
          }
          lambdaTimeTried <- 0
          lambdaTimeAccepted <- 0
        }

        # delta for value of lambda
        rate <- lambdaValueAccepted / lambdaValueTried
        if ( rate > 0.44 ) {
          deltaLambdaValue <- deltaLambdaValue * (1.0 + ((rate-0.44)/0.56) )
        } else {
          deltaLambdaValue <- deltaLambdaValue / (2.0 - rate/0.44 )
        }
        lambdaValueTried <- 0
        lambdaValueAccepted <- 0

        # delta for time of mu change events
        if ( muTimeTried > 0 ) {
          rate <- muTimeAccepted / muTimeTried
          if ( rate > 0.44 ) {
            deltaMuTime <- deltaMuTime * (1.0 + ((rate-0.44)/0.56) )
          } else {
            deltaMuTime <- deltaMuTime / (2.0 - rate/0.44 )
          }
          muTimeTried <- 0
          muTimeAccepted <- 0
        }

        # delta for values of mu
        rate <- muValueAccepted / muValueTried
        if ( rate > 0.44 ) {
          deltaMuValue <- deltaMuValue * (1.0 + ((rate-0.44)/0.56) )
        } else {
          deltaMuValue <- deltaMuValue / (2.0 - rate/0.44 )
        }
        muValueTried <- 0
        muValueAccepted <- 0

        # delta for time of massExtinction change events
        if ( massExtinctionTimeTried > 0 ) {
          rate <- massExtinctionTimeAccepted / massExtinctionTimeTried
          if ( rate > 0.44 ) {
            deltaMassExtinctionTime <- deltaMassExtinctionTime * (1.0 + ((rate-0.44)/0.56) )
          } else {
            deltaMassExtinctionTime <- deltaMassExtinctionTime / (2.0 - rate/0.44 )
          }
          massExtinctionTimeTried <- 0
          massExtinctionTimeAccepted <- 0
        }

        # delta for values of mu
        if ( massExtinctionValueTried > 0 ) {
          rate <- massExtinctionValueAccepted / massExtinctionValueTried
          if ( rate > 0.44 ) {
            deltaMassExtinctionValue <- deltaMassExtinctionValue * (1.0 + ((rate-0.44)/0.56) )
          } else {
            deltaMassExtinctionValue <- deltaMassExtinctionValue / (2.0 - rate/0.44 )
          }
          massExtinctionValueTried <- 0
          massExtinctionValueAccepted <- 0
        }
      }
    }

    # store the values
    if ( i %% THINNING == 0 & SAMPLE ) {

      sampleIndex <- length(kLambda)+1
      kLambda[sampleIndex] <- length(lambda)
      speciationRateValues[[sampleIndex]] <- lambda
      speciationRateChangeTimes[[sampleIndex]] <- lambdaChangeTimes
      kMu[sampleIndex] <- length(mu)
      extinctionRateValues[[sampleIndex]] <- mu
      extinctionRateChangeTimes[[sampleIndex]] <- muChangeTimes
      kMassExtinction[sampleIndex] <- length(pMassExtinction)
      survivalProbability[[sampleIndex]] <- pMassExtinction
      massExtinctionTime[[sampleIndex]] <- tMassExtinction
      posterior[sampleIndex] <- currentPP

      # write.table(speciationRateValues,"SpeciationRates.txt",sep="\t")
      cat(lambda,sep="\t",file="SpeciationRates.txt",append=TRUE)
      cat("\n",file="SpeciationRates.txt",append=TRUE)

      # write.table(speciationRateChangeTimes,"SpeciationRateChanges.txt",sep="\t")
      cat(lambdaChangeTimes,sep="\t",file="SpeciationRateChanges.txt",append=TRUE)
      cat("\n",file="SpeciationRateChanges.txt",append=TRUE)

      # write.table(extinctionRateValues,"ExtinctionRates.txt",sep="\t")
      cat(mu,sep="\t",file="ExtinctionRates.txt",append=TRUE)
      cat("\n",file="ExtinctionRates.txt",append=TRUE)

      # write.table(extinctionRateChangeTimes,"ExtinctionRateChanges.txt",sep="\t")
      cat(muChangeTimes,sep="\t",file="ExtinctionRateChanges.txt",append=TRUE)
      cat("\n",file="ExtinctionRateChanges.txt",append=TRUE)

      # write.table(survivalProbability,"SurvivalProbabilities.txt",sep="\t")
      cat(pMassExtinction,sep="\t",file="SurvivalProbabilities.txt",append=TRUE)
      cat("\n",file="SurvivalProbabilities.txt",append=TRUE)

      # write.table(massExtinctionTime,"MassExtinctionTimes.txt",sep="\t")
      cat(tMassExtinction,sep="\t",file="MassExtinctionTimes.txt",append=TRUE)
      cat("\n",file="MassExtinctionTimes.txt",append=TRUE)

      nSamples <- sampleIndex
      write.table(list(Iteration=nSamples*THINNING,posterior=posterior[nSamples],NumSpeciation=kLambda[nSamples],numExtinction=kMu[nSamples],numMassExtinctions=kMassExtinction[nSamples]),"samples_numCategories.txt", sep="\t", row.names = FALSE, append=TRUE, quote = FALSE, col.names = FALSE)

    }

    # Progress by convergence
    if ( i %% CONVERGENCE_FREQUENCY == 0 & SAMPLE & length(posterior) > 2 ) {

      samples <- length(posterior)

      if ( estimateNumberRateChanges == TRUE ) {
        spectralDensityLambda <- spectrum0.ar(kLambda)$spec
        spectralDensityMu     <- spectrum0.ar(kMu)$spec
        if ( spectralDensityLambda > 0 && spectralDensityMu > 0 ) {
          ess_kLambda <- (var(kLambda) * samples/spectralDensityLambda)
          ess_kMu <- (var(kMu) * samples/spectralDensityMu)
        } else {
          ess_kLambda <- ess_kMu <- Inf
        }
      } else {
        ess_kLambda <- ess_kMu <- Inf
      }

      spec <- c()
      for ( j in 1:length(speciationRateValues) ) {
        spec[j] <- speciationRateValues[[j]][1]
      }
      spectralDensity <- spectrum0.ar(spec)$spec
      if ( spectralDensity > 0 ) {
        ess_spec <- (var(spec) * samples/spectralDensity)
      } else {
        ess_spec <- Inf
      }

      ext <- c()
      for ( j in 1:length(extinctionRateValues) ) {
        ext[j] <- extinctionRateValues[[j]][1]
      }
      spectralDensity <- spectrum0.ar(ext)$spec
      if ( spectralDensity > 0 ) {
        ess_ext <- (var(ext) * samples/spectralDensity)
      } else {
        ess_ext <- Inf
      }

      if ( estimateNumberMassExtinctions == TRUE ) {
        ess_kMassExtinctions <- (var(kMassExtinction) * samples/spectrum0.ar(kMassExtinction)$spec)
      } else {
        ess_kMassExtinctions <- Inf
      }

      min.ess <- min(c(ess_kLambda,ess_kMu,ess_spec,ess_ext,ess_kMassExtinctions))

    }

    if ( SAMPLE ) {

      # Progress by iteration
      iteration.progress <- i / MAX_ITERATIONS

      # Progress by time
      time.progress <- ( (Sys.time() - startTime) / MAX_TIME )

      # Progress by ESS
      ess.progress <- min.ess / MIN_ESS

      # Max progress
      max.progress <- max(c(iteration.progress,time.progress,ess.progress))

      if ( verbose ) {
        setTxtProgressBar(bar,pmin(max.progress,1))
      }

      if( max.progress >= 1 ) {
        finished <- TRUE
      }

    } else {

      if ( verbose ){
        setTxtProgressBar(bar,pmin(1,i/burn))
      }

    }

    i <- i + 1

    if ( i == burn & !SAMPLE ) {
      startTime <- Sys.time()
      i <- 1
      SAMPLE <- TRUE
      if ( verbose ) {
        cat("\n\nRunning the chain ... \n")
        cat("0--------25--------50--------75--------100\n")
        bar <- txtProgressBar(style=1,width=42)
      }
    }

  }

  cat("\n")
  setwd(file.path(orgDir))

}


globalBiDe.analysis <- tess.analysis













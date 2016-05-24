################################################################################
#
# tess.likelihood.rateshift.R
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
# @brief Computation of the likelihood for a given tree under a birth-death rate-shift model (i.e. piecewise constant rates).
#
# @date Last modified: 2015-05-28
# @author Sebastian Hoehna
# @version 2.0
# @since 2012-09-22, version 1.3
#
# @param    times                                         vector        vector of branching times
# @param    times                                         vector        branching times
# @param    lambda                                        vector        speciation rates
# @param    mu                                            vector        extinction rates
# @param    rateChangeTimesLambda                         vector        speciation rates
# @param    rateChangeTimesMu                             vector        extinction rates
# @param    massExtinctionTimes                           vector        time at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified|age
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @param    log                                           boolean       likelhood in log-scale?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

tess.likelihood.rateshift <- function( times,
                                             lambda,
                                             mu,
                                             rateChangeTimesLambda = c(),
                                             rateChangeTimesMu = c(),
                                             massExtinctionTimes = c(),
                                             massExtinctionSurvivalProbabilities = c(),
                                             missingSpecies = c(),
                                             timesMissingSpecies = c(),
                                             samplingStrategy = "uniform",
                                             samplingProbability = 1.0,
                                             MRCA=TRUE,
                                             CONDITION="survival",
                                             log=TRUE) {
  
  if ( length(lambda) != (length(rateChangeTimesLambda)+1) || length(mu) != (length(rateChangeTimesMu)+1) ) {
    stop("Number of rate-change times needs to be one less than the number of rates!")
  }

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equal the number of mass-extinction survival probabilities!")
  }

  if ( length(missingSpecies) != length(timesMissingSpecies) ) {
    stop("Vector holding the missing species must be of the same size as the intervals when the missing speciation events happend!")
  }

  if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "taxa" ) {
    stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|taxa.")
  }

  if ( samplingStrategy != "uniform" && samplingStrategy != "diversified") {
    stop("Wrong choice of argument for \"samplingStrategy\". Possible option are uniform|diversified.")
  }

    # make sure the times and values are sorted
  if ( length(rateChangeTimesLambda) > 0 ) {
    sortedRateChangeTimesLambda <- sort( rateChangeTimesLambda )
    lambda <- c(lambda[1], lambda[ match(sortedRateChangeTimesLambda,rateChangeTimesLambda)+1 ] )
    rateChangeTimesLambda <- sortedRateChangeTimesLambda
  }
  if ( length(rateChangeTimesMu) > 0 ) {
    sortedRateChangeTimesMu <- sort( rateChangeTimesMu )
    mu <- c(mu[1], mu[ match(sortedRateChangeTimesMu,rateChangeTimesMu)+1 ] )
    rateChangeTimesMu <- sortedRateChangeTimesMu
  }
  if ( length(massExtinctionTimes) > 0 ) {
    sortedMassExtinctionTimes <- sort( massExtinctionTimes )
    massExtinctionSurvivalProbabilities <- massExtinctionSurvivalProbabilities[ match(sortedMassExtinctionTimes,massExtinctionTimes) ]
    massExtinctionTimes <- sortedMassExtinctionTimes
  }
  
  # join the times of the rate changes and the mass-extinction events
  if ( length( rateChangeTimesLambda ) > 0 ||  length( rateChangeTimesMu ) > 0 || length( massExtinctionTimes ) > 0 ) {
    changeTimes <- sort( unique( c( rateChangeTimesLambda, rateChangeTimesMu, massExtinctionTimes ) ) )
  } else {
    changeTimes <- c()
  }
  speciation <- rep(NaN,length(changeTimes)+1)
  extinction <- rep(NaN,length(changeTimes)+1)
  mep <- rep(NaN,length(changeTimes))
  speciation[1] <- lambda[1]
  if ( length(lambda) > 1 ) {
    speciation[ match(rateChangeTimesLambda,changeTimes)+1 ] <- lambda[ 2:length(lambda) ]
  }
  extinction[1] <- mu[1]
  if ( length(mu) > 1 ) {
    extinction[ match(rateChangeTimesMu,changeTimes)+1 ] <- mu[ 2:length(mu) ]
  }
  if ( length( massExtinctionSurvivalProbabilities ) > 0 ) {
    mep[ match(massExtinctionTimes,changeTimes) ] <- massExtinctionSurvivalProbabilities[ 1:length(massExtinctionSurvivalProbabilities) ]
  }
  for ( i in seq_len(length(changeTimes)) ) {
    if ( is.null(speciation[i+1]) || !is.finite(speciation[i+1]) ) {
      speciation[i+1] <- speciation[i]
    }
    if ( is.null(extinction[i+1]) || !is.finite(extinction[i+1]) ) {
      extinction[i+1] <- extinction[i]
    }
    if ( is.null(mep[i]) || !is.finite(mep[i]) ) {
      mep[i] <- 1.0
    }
  }

  
  rateChangeTimes <- changeTimes
  massExtinctionTimes <- changeTimes

  lambda <- speciation
  mu <- extinction
  massExtinctionSurvivalProbabilities <- mep

  
  PRESENT <- max(times)
  nTaxa <- length(times) + 1
  
  times <- PRESENT - sort(times,decreasing=TRUE)

  # if we condition on the MRCA, then we need to remove the root speciation event
  if ( MRCA == TRUE ) {
    times <- times[-1]
  }

  # set the uniform taxon sampling probability
  if (samplingStrategy == "uniform") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }
  
  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,PRESENT,log=TRUE)
  
  # multiply the probability of a descendant of the initial species
  lnl <- lnl + tess.equations.p1.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl
  } 

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl + tess.equations.pN.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,rho,nTaxa,0,PRESENT,SURVIVAL=TRUE,MRCA,log=TRUE)

  # if we assume diversified sampling, we need to multiply with the probability that all missing species happened after the last speciation event
  if ( samplingStrategy == "diversified" ) {
    # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
    lastEvent <- times[length(times)]
    p_0_T <- 1.0 - tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT,log=FALSE) * exp((mu-lambda)*PRESENT)
    p_0_t <- 1.0 - tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT,log=FALSE)*exp((mu-lambda)*(PRESENT-lastEvent))
    F_t <- p_0_t / p_0_T
    # get an estimate of the actual number of taxa
    m <- round(nTaxa / samplingProbability)
    # remove the number of species that we started with
    k <- 1
    if ( MRCA == TRUE ) k <- 2
    lnl <- lnl + (m-nTaxa) * log(F_t) + lchoose(m-k,nTaxa-k)
  }

  # multiply the probability for the missing species
  if ( length(missingSpecies) > 0 ) {

    # compute the rate
    prev_time <- 0
    rate <- 0
    # add mass-extinction
    for (j in seq_len(length(rateChangeTimes)) ) {
      rate <- rate + ifelse( PRESENT >= rateChangeTimes[j], (mu[j] - lambda[j])*(rateChangeTimes[j]-prev_time) - log(massExtinctionSurvivalProbabilities[j]), 0 )
      prev_time <- ifelse( PRESENT >= rateChangeTimes[j], rateChangeTimes[j], 0) 
    }
    # add the final rate interval
    rate <- rate + ifelse( PRESENT > prev_time, (mu[length(mu)] - lambda[length(lambda)])*(PRESENT-prev_time), 0 )

    # add sampling
    rate <- rate - log(samplingProbability)
    
    p_0_T <- 1.0 - exp( tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT,log=TRUE) + rate )

    # now iterate over the vector of missing species per interval
    lastEvent <- timesMissingSpecies
    
    # compute the rate
    prev_time <- lastEvent
    rate <- 0
    # add mass-extinction
    for (j in seq_len(length(rateChangeTimes)) ) {
      rate <- rate + ifelse( lastEvent < rateChangeTimes[j] & PRESENT >= rateChangeTimes[j], (mu[j] - lambda[j])*(rateChangeTimes[j]-prev_time) - log(massExtinctionSurvivalProbabilities[j]), 0 )
      prev_time <- ifelse( lastEvent < rateChangeTimes[j] & PRESENT >= rateChangeTimes[j], rateChangeTimes[j], lastEvent) 
    }
    # add the final rate interval
    rate <- rate + ifelse( PRESENT > prev_time, (mu[length(mu)] - lambda[length(lambda)])*(PRESENT-prev_time), 0 )

    # add sampling
    rate <- rate - log(samplingProbability)

  
    p_0_t <- 1.0 - exp( tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT,log=TRUE) + rate )
    log_F_t <- log(p_0_t) - log(p_0_T)
    # get an estimate of the actual number of taxa
    m <- missingSpecies
    # remove the number of species that we started with

    lnl <- lnl + sum( m * log_F_t ) #+ lchoose(m-k,nTaxa-k)
  
  }

  if ( length(rateChangeTimes) > 0 ) {
    speciation <- function(times) {
      idx <- findInterval(times,rateChangeTimes)+1
      idx[ idx > length(lambda) ] <- length(lambda)
      return ( lambda[idx] )
    }
  } else {
    speciation <- function(times) rep(lambda[1],length(times))
  }
    
  # multiply the probability for each speciation time  
  lnl <- lnl + sum( log(speciation(times) ) ) + sum(tess.equations.p1.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,rho,times,PRESENT,log=TRUE))
  
  if (is.nan(lnl)) lnl <- -Inf

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }

  return (lnl)
}






################################################################################
#
# tess.likelihood.R
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
# @brief Computation of the likelihood for a given tree.
#
# @date Last modified: 2015-05-28
# @author Sebastian Hoehna
# @version 1.2
# @since 2012-09-22, version 1.0
#
# @param    times                                         vector        vector of branching times
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified|age
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @param    log                                           boolean       likelhood in log-scale?

# @return                                                 scalar        probability of the speciation times
#
################################################################################

tess.likelihood <- function(times,
                            lambda,
                            mu,
                            massExtinctionTimes=c(),
                            massExtinctionSurvivalProbabilities=c(),
                            missingSpecies = c(),
                            timesMissingSpecies = c(),
                            samplingProbability=1.0,
                            samplingStrategy="uniform",
                            MRCA=TRUE,
                            CONDITION="survival",
                            log=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( CONDITION != "time" && CONDITION != "survival" && CONDITION != "taxa" ) {
    stop("Wrong choice of argument for \"CONDITION\". Possible option are time|survival|taxa.")
  }

  if ( samplingStrategy != "uniform" && samplingStrategy != "diversified" ) {
    stop("Wrong choice of argument for \"samplingStrategy\". Possible option are uniform|diversified.")
  }

  if ( (!is.numeric(lambda) && !inherits(lambda, "function")) || (!is.numeric(mu) && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }

  # test if we got constant values for the speciation and extinction rates
  if ( is.numeric(lambda) && is.numeric(mu) ) {
    # call computations for constant rates (much faster)
    p <- tess.likelihood.constant(times,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,missingSpecies,timesMissingSpecies,samplingProbability,samplingStrategy,MRCA,CONDITION,log)
    return (p)
  } else  {

    # convert the speciation rate into a function if necessary
    if ( is.numeric(lambda) ) {
      speciation <- function (x) rep(lambda,length(x))
    } else {
      speciation <- lambda
    }
    # convert the extinction rate into a function if necessary
    if ( is.numeric(mu) ) {
      extinction <- function (x) rep(mu,length(x))
    } else {
      extinction <- mu
    }

    p <- tess.likelihood.function(times,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,missingSpecies,timesMissingSpecies,c(),samplingProbability,samplingStrategy,MRCA,CONDITION,log)

    return (p)
  }

}



################################################################################
#
# @brief Computation of the likelihood for a given tree.
#
# Here we use equation (6) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
# For the diversified taxon sampling, see Equation (6) Hoehna, S., 2013, Likelihood Inference of non-constant Diversification Rates with Incomplete Taxon Sampling
# and equation (5) of Hoehna et al. (2011) "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
#
# @date Last modified: 2015-05-28
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-22, version 1.0
#
# @param    times                                         vector        vector of branching times
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

tess.likelihood.constant <- function(times,
                                     lambda,
                                     mu,
                                     massExtinctionTimes,
                                     massExtinctionSurvivalProbabilities,
                                     missingSpecies = c(),
                                     timesMissingSpecies = c(),
                                     samplingProbability,
                                     samplingStrategy,
                                     MRCA,
                                     CONDITION,
                                     log) {

  # check for sensible parameter values
  if ( lambda <= 0 || mu < 0 || samplingProbability <= 0 || samplingProbability > 1.0) {
    stop("Invalid parameter values for lambda and mu!")
  }

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
  if ( CONDITION == "survival" )    lnl <- - tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,PRESENT,log=TRUE)

  # multiply the probability of a descendant of the initial species
  lnl <- lnl + tess.equations.p1.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl
  }

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl - tess.equations.pN.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,0,PRESENT,SURVIVAL=FALSE,MRCA,log=TRUE)

  # if we assume diversified sampling, we need to multiply with the probability that all missing species happened after the last speciation event
  if ( samplingStrategy == "diversified" ) {
    # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"
    lastEvent <- times[length(times)]
    p_0_T <- 1.0 - tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT,log=FALSE) * exp((mu-lambda)*PRESENT)
    p_0_t <- 1.0 - tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT,log=FALSE) * exp((mu-lambda)*(PRESENT-lastEvent))
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
    # We use equation (5) of Hoehna et al. "Inferring Speciation and Extinction Rates under Different Sampling Schemes"

    # now iterate over the vector of missing species per interval
    lastEvent <- timesMissingSpecies

    p_0_T <- 1.0 - exp( tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,0,PRESENT,PRESENT,log=TRUE) + ((mu-lambda)*PRESENT) )
    p_0_t <- 1.0 - exp( tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,1.0,lastEvent,PRESENT,PRESENT,log=TRUE) + ((mu-lambda)*(PRESENT-lastEvent)) )

    log_F_t <- log(p_0_t) - log(p_0_T)
    # get an estimate of the actual number of taxa
    m <- missingSpecies
    # remove the number of species that we started with

    lnl <- lnl + sum( m * log_F_t ) #+ lchoose(m-k,nTaxa-k)

  }


  # multiply the probability for each speciation time
  if ( length(times) > 0 ) {
    lnl <- lnl + length(times)*log(lambda) + sum(tess.equations.p1.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,rho,times,PRESENT,log=TRUE))
  }

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }

  if (is.nan(lnl)) lnl <- -Inf

  return (lnl)
}



################################################################################
#
# @brief Computation of the likelihood for a given tree.
#
# Here we use equation (6) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
# For the diversified taxon sampling, see Equation (6) in Hoehna, S., 2014, Likelihood Inference of non-constant Diversification Rates with Incomplete Taxon Sampling
#
# @date Last modified: 2015-05-28
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-22, version 1.0
#
# @param    times                                         vector        vector of branching times
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        Which strategy was used to obtain the samples (taxa). Options are: uniform|diversified
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    CONDITITON                                    string        do we condition the process on nothing|survival|taxa?
# @return                                                 scalar        probability of the speciation times
#
################################################################################

tess.likelihood.function <- function(times,
                                     lambda,
                                     mu,
                                     massExtinctionTimes,
                                     massExtinctionSurvivalProbabilities,
                                     missingSpecies = c(),
                                     timesMissingSpecies = c(),
                                     t.crit,
                                     samplingProbability,
                                     samplingStrategy,
                                     MRCA,
                                     CONDITION,
                                     log) {

  lnl <- -Inf

#   times <- as.numeric(branching.times(tree))
  PRESENT <- max(times)
  nTaxa <- length(times) + 1

  times <- PRESENT - sort(times,decreasing=TRUE)

  # if we condition on the MRCA, then we need to remove the root speciation event
  if (MRCA == TRUE) {
    times <- times[-1]
  }

  # set the uniform taxon sampling probability
  if (samplingStrategy == "uniform") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }

  # prepare the integrals
  tryCatch({
    approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,PRESENT,c(t.crit,times),TRUE)

  # initialize the log likelihood
  lnl <- 0

  # what do we condition on?
  # did we condition on survival?
  if ( CONDITION == "survival" || CONDITION == "taxa" )    lnl <- - tess.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,rho,0,PRESENT,PRESENT,log=TRUE)


  # multiply the probability of a descendant of the initial species
  lnl <- lnl + tess.equations.p1.fastApprox(approxFuncs$r,approxFuncs$s,rho,0,PRESENT,log=TRUE)

  # add the survival of a second species if we condition on the MRCA
  if ( MRCA == TRUE ) {
    lnl <- 2*lnl
  }

  # did we condition on observing n species today
  if ( CONDITION == "taxa" )    lnl <- lnl - tess.equations.pN.fastApprox(approxFuncs$r,approxFuncs$s,rho,nTaxa,0,PRESENT,SURVIVAL=FALSE,MRCA,log=TRUE)

  # if we assume diversified sampling, we need to multiply with the probability that all missing species happened after the last speciation event
  if ( samplingStrategy == "diversified" ) {
    # We use equation (6) of Hoehna (2014). " Likelihood Inference of non-constant Diversification Rates with Incomplete Taxon Sampling"
    lastEvent <- times[length(times)]

    p_0_T <- 1.0 - tess.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,1.0,0,PRESENT,PRESENT,log=FALSE) * exp(approxFuncs$r(PRESENT))
    p_0_t <- (1.0 - tess.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,1.0,lastEvent,PRESENT,PRESENT,log=FALSE) * exp(approxFuncs$r(PRESENT) - approxFuncs$r(lastEvent)))
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
    # We use equation (6) of Hoehna (2014). " Likelihood Inference of non-constant Diversification Rates with Incomplete Taxon Sampling"

    # now iterate over the vector of missing species per interval
    lastEvent <- timesMissingSpecies

    p_0_T <- 1.0 - exp( tess.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,1.0,0,PRESENT,PRESENT,log=TRUE) + (approxFuncs$r(PRESENT)) )
    p_0_t <- 1.0 - exp( tess.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,1.0,lastEvent,PRESENT,PRESENT,log=TRUE) + (approxFuncs$r(PRESENT) - approxFuncs$r(lastEvent)) )
    log_F_t <- log(p_0_t) - log(p_0_T)
    # get an estimate of the actual number of taxa
    m <- missingSpecies
    # remove the number of species that we started with

    lnl <- lnl + sum( m * log_F_t ) #+ lchoose(m-k,nTaxa-k)

  }


  # multiply the probability for each speciation time
  if ( length(times) > 0 ) {
    lnl <- lnl + sum(log(lambda(times)) + tess.equations.p1.fastApprox(approxFuncs$r,approxFuncs$s,rho,times,PRESENT,log=TRUE))
  }

  }, warning = function(w) {  }, error = function(e) { })

  if ( log == FALSE ) {
    lnl <- exp(lnl)
  }

  if (is.nan(lnl)) lnl <- -Inf

  return (lnl)


}





################################################################################
#
# tess.sim.age.R
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
# @brief Simulate a tree for a given age (either using numerical integration
#        or the analytical solution for constant rates).
#
# @date Last modified: 2014-10-05
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-23, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        the incomplete taxon sampling strategy, either uniform or diversified
# @param    maxTaxa                                       scalar        the maximum number of taxa that will be simulated. Once the limit has been reached it will terminate.
# @param    MRCA                                          boolean       start with a single or two species?
# @return                                                 phylo         a simulated tree
#
################################################################################
tess.sim.age <- function(n,age,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,samplingStrategy="uniform",maxTaxa=Inf,MRCA=TRUE) {

  if ( length(massExtinctionTimes) != length(massExtinctionSurvivalProbabilities) ) {
    stop("Number of mass-extinction times needs to equals the number of mass-extinction survival probabilities!")
  }

  if ( samplingStrategy != "uniform" && samplingStrategy != "diversified" ) {
    stop("Wrong choice of argument for \"samplingStrategy\". Possible option are uniform|diversified.")
  }

  if ( (!is.numeric(lambda) && !inherits(lambda, "function")) || (!is.numeric(mu) && !inherits(mu, "function"))) {
    stop("Unexpected parameter types for lambda and mu!")
  }
  
  # test if we got constant values for the speciation and extinction rates
  if ( is.numeric(lambda) && is.numeric(mu) ) {
    # call simulation for constant rates (much faster)
    trees <- tess.sim.age.constant(n,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,maxTaxa,MRCA)
    return (trees)
  } else {
    
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

    trees <- tess.sim.age.function(n,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,maxTaxa,MRCA)
    return (trees)
  }

}
 


################################################################################
# 
# @brief Simulate a tree for a given age.
#
# 1) Draw n times the number of taxa at the present time (Equation (10) in Hoehna,
# Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes.
# 2013, Bioinformatics, 29:1367-1374 ).
# 2) For each nTaxa, draw simulate a tree using the function tess.sim.taxa.age.constant
#
# @date Last modified: 2014-01-14
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        the incomplete taxon sampling strategy, either uniform or diversified
# @param    maxTaxa                                       scalar        the maximum number of taxa that will be simulated. Once the limit has been reached it will terminate.
# @param    MRCA                                          boolean       start with a single or two species?
# @return                                                 phylo         a simulated tree
#
################################################################################
tess.sim.age.constant <- function(n,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,maxTaxa,MRCA) {

  # check for sensible parameter values
  if ( lambda <= 0 || mu < 0 || samplingProbability <= 0 || samplingProbability > 1.0) {
    stop("Invalid parameter values for lambda and mu!")
  }

  
  # precompute the probabilities
  p_s <- tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,age,age,log=FALSE)
  r   <- (mu-lambda)*age - log(samplingProbability)
  for (j in seq_len(length(massExtinctionTimes)) ) {
    cond <-  (0 < massExtinctionTimes[j]) & (age >= massExtinctionTimes[j])
    r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
  }

  # randomly draw a number of taxa
  m <- rgeom(n, p_s * exp(r)) + 1
  if ( MRCA == TRUE ) m <- m + rgeom(n, p_s * exp(r)) + 1
  
  # for safety we set default values
  m[is.na(m)] <- maxTaxa

  m <- pmin(maxTaxa,m)
  
  trees <- list()
  # for each simulation
  for ( i in 1:n ) {
    nTaxa <- m[i]

    if ( !is.finite( nTaxa ) ) stop( sprintf("Simulation for age = %.2f, lambda = %.2f and mu = %.2f resulted in a non-finite number of taxa.", age, lambda, mu) )

    # check if we actually have a tree
    if ( nTaxa < 2 ) {
      tree <- 1
    } else {

      # delegate the call to the simulation condition on both, age and nTaxa
      tree <- tess.sim.taxa.age.constant(1,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)[[1]]
    }

    trees[[i]] <- tree
  } # end for each simulation

  return (trees)
}



################################################################################
# 
# @brief Simulate a tree for a given age.
#
# 1) Draw n times the number of taxa at the present time (Equation (10) in Hoehna,
# Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes.
# 2013, Bioinformatics, 29:1367-1374 ).
# 2) For each nTaxa, draw simulate a tree using the function tess.sim.taxa.age.constant
#
# @date Last modified: 2014-01-14
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        the incomplete taxon sampling strategy, either uniform or diversified
# @param    maxTaxa                                       scalar        the maximum number of taxa that will be simulated. Once the limit has been reached it will terminate.
# @param    MRCA                                          boolean       start with a single or two species?
# @return                                                 phylo         a simulated tree
#
################################################################################
tess.sim.age.rateshift <- function( n,
                                    age,
                                    lambda,
                                    mu,
                                    rateChangeTimesLambda = c(),
                                    rateChangeTimesMu = c(),
                                    massExtinctionTimes = c(),
                                    massExtinctionSurvivalProbabilities = c(),
                                    samplingProbability,
                                    samplingStrategy,
                                    maxTaxa,
                                    MRCA) {

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

  
  # precompute the probabilities
  p_s <- tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,0,age,age,log=FALSE)
  r   <- (mu-lambda)*age - log(samplingProbability)
  for (j in seq_len(length(massExtinctionTimes)) ) {
    cond <-  (0 < massExtinctionTimes[j]) & (age >= massExtinctionTimes[j])
    r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
  }

  # randomly draw a number of taxa
  m <- rgeom(n, p_s * exp(r)) + 1
  if ( MRCA == TRUE ) m <- m + rgeom(n, p_s * exp(r)) + 1

  m[is.na(m)] <- maxTaxa
  m <- pmin(maxTaxa,m)
  
  trees <- list()
  # for each simulation
  for ( i in 1:n ) {
    nTaxa <- m[i]

    if ( !is.finite( nTaxa ) ) stop( sprintf("Simulation for age = %.2f, lambda = %.2f and mu = %.2f resulted in a non-finite number of taxa.", age, lambda, mu) )

    # check if we actually have a tree
    if ( nTaxa < 2 ) {
      tree <- 1
    } else {

      # delegate the call to the simulation condition on both, age and nTaxa
      tree <- tess.sim.taxa.age.rateshift(1,nTaxa,age,lambda,mu,rateChangeTimes,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)[[1]]
    }

    trees[[i]] <- tree
  } # end for each simulation

  return (trees)
}




################################################################################
# 
# @brief Simulate a tree for a given age. This is the fast approximation
#        procedure using precomputed integral functions.
#
# 1) Draw n times the number of taxa at the present time (Equation (10) in Hoehna,
# Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes.
# 2013, Bioinformatics, 29:1367-1374 ).
# 2) For each nTaxa, draw simulate a tree using the function tess.sim.taxa.age.function
#
# @date Last modified: 2014-10-05
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    age                                           scalar        time of the process
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    samplingStrategy                              string        the incomplete taxon sampling strategy, either uniform or diversified
# @param    maxTaxa                                       scalar        the maximum number of taxa that will be simulated. Once the limit has been reached it will terminate.
# @param    MRCA                                          boolean       start with a single or two species?
# @return                                                 phylo         a simulated tree
#
################################################################################
tess.sim.age.function <- function(n,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,maxTaxa,MRCA) {

  # set the uniform taxon sampling probability
  if (samplingStrategy == "uniform") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }

  # approximate the rate integral and the survival probability integral for fas computations
  approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,age,c())

  # precompute the probabilities
  p_s <- tess.equations.pSurvival.fastApprox(approxFuncs$r,approxFuncs$s,rho,0,age,age,log=FALSE)
  r   <- approxFuncs$r(age) - log(rho)

  # randomly draw a number of taxa
  nTaxa <- rgeom(n, p_s * exp(r)) + 1
  if ( MRCA == TRUE ) nTaxa <- nTaxa + rgeom(n, p_s * exp(r)) + 1
  
  nTaxa[is.na(nTaxa)] <- maxTaxa
  nTaxa <- pmin(maxTaxa,nTaxa)

  trees <- list()
  # for each simulation
  for ( i in 1:n ) {

    # check if we actually have a tree
    if ( nTaxa[i] < 2 ) {
      tree <- 1
    } else {

      # delegate the call to the simulation condition on both, age and nTaxa
      tree <- tess.sim.taxa.age.function(1,nTaxa[i],age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA=MRCA,approxFuncs$r,approxFuncs$s)[[1]]
    }

    trees[[i]] <- tree
  } # end for each simulation
  
  return (trees)
    
}




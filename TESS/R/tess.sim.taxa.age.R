################################################################################
#
# globalBiDe.equations.R
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
# @brief Simulate a tree for a given number of taxa.
#
# @date Last modified: 2012-11-05
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
tess.sim.taxa.age <- function(n,nTaxa,age,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,samplingStrategy="uniform",MRCA=TRUE) {

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
    trees <- tess.sim.taxa.age.constant(n,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)
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

    approxFuncs <- tess.prepare.pdf(speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,age,c())
    
    trees <- tess.sim.taxa.age.function(n,nTaxa,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,approxFuncs$r,approxFuncs$s)
    
    return (trees)
  } 

}


################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree.
#
#
# Simulate the n speciation times using the inverse cdf given in Equation (9)
# in Hoehna, Fast simulation of reconstructed phylogenies under global,
# time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 .
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @return                                                 phylo         a random tree
#
################################################################################
tess.sim.taxa.age.constant <- function(n,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA) {

  # check for sensible parameter values
  if ( lambda <= 0 || mu < 0 || samplingProbability <= 0 || samplingProbability > 1.0) {
    stop("Invalid parameter values for lambda and mu!")
  }

  # if we have mass-extinction events, then we use the function-integration approach
  # because the closed form solutions only apply to models without mass-extinctions.
  if ( length(massExtinctionTimes) > 0 ) {
    speciation <- function(x) lambda
    extinction <- function(x) mu
    approxFuncs <- tess.prepare.pdf(speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,age,c())    
    trees <- tess.sim.taxa.age.function(n,nTaxa,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,approxFuncs$r,approxFuncs$s)
    return (trees)
  } else {

    # set the uniform taxon sampling probability
    if (samplingStrategy == "uniform") {
      rho <- samplingProbability
    } else {
      rho <- 1.0
    }
    
    # if we use diversified sampling then we just sample m instead of n speciation events
    if (samplingStrategy == "diversified") {
      m <- nTaxa
      nTaxa <- round(nTaxa / samplingProbability)
      nTaxa <- pmin(nTaxa,10000)
    }
    
#    b <- rho * lambda
#    d <- mu - lambda *(1-rho)
    b <- lambda
    d <- mu

    # start to simulate n trees
    trees <- list()
    for ( j in 1:n ) {

      # now draw the random time speciation times
      times         <- c()
      x <- 1
      if (MRCA == TRUE) x <- 2
      if ( (nTaxa-x) > 0) {
        
        # for each speciation time
        u           <- runif(nTaxa-x,0,1)
#        times       <- 1/(b-d)*log((b - d*exp((-b+d)*age) -d*(1-exp((-b+d)*age)) *u )/(b - d*exp((-b+d)*age) -b*(1.0-exp((-b+d)*age)) *u )   )  
#        times       <- log( 1.0 - u*(1-exp((-b+d)*age)) ) / (d-b)
        times       <- age - ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(rho*b+(b*(1-rho)-d)*exp((d-b)*age) ) ) ) - (b*(1-rho)-d) ) / (rho * b) ) + (d-b)*age )  /  (d-b) 
# now we have the vector of speciation times, which just need to be converted to a phylogeny
        times         <- c(sort(times,decreasing=FALSE),age)
        
        # if we use diversified sampling we need to remove the last n-m species
        if (samplingStrategy == "diversified") {
          if ( m < nTaxa ) {
            times <- times[-(1:(nTaxa-m))]
          }
        }
      }
      
      if ( length(times) < 1 ) {
        times         <- c(rep(age/2,2-x),age)
      }
      
      tree          <- tess.create.phylo(times, root = !MRCA)

      trees[[j]]    <- tree
    } # end for each simulation

    return (trees)
  }
}


################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree.
#
#
# Simulate the n speciation times using the inverse cdf given in Equation (9)
# in Hoehna, Fast simulation of reconstructed phylogenies under global,
# time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 .
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @return                                                 phylo         a random tree
#
################################################################################
tess.sim.taxa.age.rateshift <- function(n,nTaxa,age,lambda,mu,rateChangeTimes,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA) {

  # if we have mass-extinction events, then we use the function-integration approach
  # because the closed form solutions only apply to models without mass-extinctions.
 
  speciation <- function(x) lambda[findInterval(t,rateChangeTimes)]
  extinction <- function(x) mu[findInterval(t,rateChangeTimes)]
  approxFuncs <- tess.prepare.pdf(speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,age,rateChangeTimes)    
  trees <- tess.sim.taxa.age.function(n,nTaxa,age,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,approxFuncs$r,approxFuncs$s)
  return (trees)
 
}




################################################################################
# 
# @brief Simulate a tree conditioned on the number of taxa and
#        the age of the tree. This is the fast approximation
#        procedure using precomputed integral functions.
#
# Simulate the n speciation times using the inverse cdf given in Equation (9)
# in Hoehna, Fast simulation of reconstructed phylogenies under global,
# time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 .
#
# @date Last modified: 2012-12-17
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    age                                           scalar        the age
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    approx                                        boolean       should linear function approximation be used
# @return                                                 phylo         a random tree
#
################################################################################
tess.sim.taxa.age.function <- function(n,nTaxa,age,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,r,s) {
  
  # set the uniform taxon sampling probability
  if (samplingStrategy == "uniform") {
    rho <- samplingProbability
  } else {
    rho <- 1.0
  }


  n2 <- 1001
  times <- seq(0, age, length=n2)
  const_zz <- 1.0 - tess.equations.pSurvival.fastApprox(r,s,rho,0,age,age,log=FALSE) * exp(r(age) - log(rho))
  zz <- 1.0 - ( (1.0 - tess.equations.pSurvival.fastApprox(r,s,rho,times,age,age,log=FALSE) * exp(r(age) - log(rho) - r(times))) / const_zz )
  zz[n2] <- 1.0
  icdf <- approxfun(zz, times) ## Interpolate
  speciationEventSim <- function(n)  icdf(runif(n))
  
  # if we use diversified sampling then we just sample m instead of n speciation events
  if (samplingStrategy == "diversified") {
    m <- nTaxa
    nTaxa <- round(nTaxa / samplingProbability)
  }

  # start to simulate n trees
  trees <- list()
  for ( j in 1:n ) {
    
    # now draw the random time speciation times
    times         <- c()
    x <- 1
    if (MRCA == TRUE) x <- 2
    if ( (nTaxa-x) > 0) {
      times <- speciationEventSim(nTaxa-x)

      # now we have the vector of speciation times, which just need to be converted to a phylogeny
      times         <- c(age - sort(times,decreasing=TRUE),age)
      
      # if we use diversified sampling we need to remove the last n-m species
      if (samplingStrategy == "diversified") {
        if ( m < nTaxa ) {
          times <- times[-(1:(nTaxa-m))]
        }
      }
    }
    
    if ( length(times) < 1 ) {
        times         <- c(rep(age/2,2-x),age)
    }

    # create the tree
    tree          <- tess.create.phylo( times, root = !MRCA )

    trees[[j]]    <- tree
  } #end for each simulation
  
  return (trees)
}

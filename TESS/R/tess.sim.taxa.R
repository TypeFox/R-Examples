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
# @brief Simulate a tree for a given number of taxa.
#
# @date Last modified: 2012-11-05
# @author Sebastian Hoehna
# @version 1.0
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    max                                           scalar        the maximal time for the age
# @param    lambda                                        function      speciation rate function
# @param    mu                                            function      extinction rate function
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    MRCA                                          boolean       do we start the tree with the MRCA (two species)?
# @param    t_crit                                        vector        critical times when jumps in the rate functions occur
# @return                                                 list          list of random trees (type phylo)
#
################################################################################
tess.sim.taxa <- function(n,nTaxa,max,lambda,mu,massExtinctionTimes=c(),massExtinctionSurvivalProbabilities=c(),samplingProbability=1.0,samplingStrategy="uniform",SURVIVAL=TRUE,MRCA=TRUE,t_crit=c()) {

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
    trees <- tess.sim.taxa.constant(n,nTaxa,max,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,SURVIVAL,MRCA)
    return (trees)
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
    
    trees <- tess.sim.taxa.function(n,nTaxa,max,speciation,extinction,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,SURVIVAL,MRCA,t_crit)
    return (trees)
  }

}



################################################################################
# 
# @brief Simulate a tree for a given number of taxa.
#
# 1) Simulate the time of the process using Monte Carlo sampling, see Equation (11)
# in Hoehna, Fast simulation of reconstructed phylogenies under global,
# time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 .
# 2) Simulate the tree by calling tess.sim.taxa.age.constant
# Note: The sampling strategy does not affect the probability of the age of the tree
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    n                                             scalar        number of simulations
# @param    nTaxa                                         scalar        number of taxa at present
# @param    max                                           scalar        the maximal time for the age
# @param    lambda                                        scalar        speciation rate function
# @param    mu                                            scalar        extinction rate function
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @return                                                 list          list of random trees (type phylo)
#
################################################################################
tess.sim.taxa.constant <- function(n,nTaxa,max,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,SURVIVAL,MRCA) {

  # check for sensible parameter values
  if ( lambda <= 0 || mu < 0 || samplingProbability <= 0 || samplingProbability > 1.0) {
    stop("Invalid parameter values for lambda and mu!")
  }

# We will always use the numerical procedures to sample form the root age
#  if ( length(massExtinctionTimes) >= 0 ) {
    # compute the cumulative distribution function for the time of the process
    pdf <- function(x) tess.equations.pN.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,nTaxa,0,x,SURVIVAL,MRCA,log=FALSE)

    # use the inverse-cdf to sample
    repeat { # find a better interval
      obj.pdf <- function(t, state, pars) list(pdf(t))
      ## This step is slow for large n - perhaps 1/2s for 1000 points
      n1 <- 11
      times <- seq(0, max, length=n1)
      zz <- lsoda(0, times, obj.pdf, tcrit=max)[,2]
      m <- min(which( (zz[n1] - zz)/zz[n1] < 1E-5))
      if ( m > (n1/2) ) {
        # now we do it properly
        n2 <- 101
        times <- seq(0, max, length=n2)
        zz <- lsoda(0, times, obj.pdf, tcrit=max)[,2]
        zz <- zz / zz[n2]  # normalize
        break
      } else {
        max <- max/2.0
      }
    }

    icdf <- approxfun(zz, times) ## Interpolate
    simT <- function(n)  icdf(runif(n))
    T <- simT(n)
#  } else {
#    m <- rep(nTaxa / samplingProbability, n)
#    u <- runif(n,0,1)
#    T <- log((-lambda - lambda * u^(1/m) + mu * u^(1/m) + lambda * u^(1/m))/(lambda * (-1 + u^(1/m)))) / (lambda - mu)
#  }
    
  trees <- list()
  # for each simulation
  for ( i in 1:n ) {

    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- tess.sim.taxa.age.constant(1,nTaxa,T[i],lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA)
  
    # add the new tree to the list
    trees[[i]] <- tree[[1]]
  }

  return (trees)
}



################################################################################
# 
# @brief Simulate a tree for a given number of taxa.
#
# 1) Simulate the time of the process using Monte Carlo sampling, see Equation (11)
# in Hoehna, Fast simulation of reconstructed phylogenies under global,
# time-dependent birth-death processes. 2013, Bioinformatics, 29:1367-1374 .
# 2) Simulate the tree by calling sim.globalBiDe.taxa.age.function
# Note: The sampling strategy does not affect the probability of the age of the tree
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda        function      speciation rate function
# @param    mu            function      extinction rate function
# @param    nTaxa         scalar        number of taxa in the tree at the present time
# @return                 phylo         a random tree
#
################################################################################
tess.sim.taxa.function <- function(n,nTaxa,max,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,SURVIVAL,MRCA,t_crit=c()) {

  repeat { # find a better interval
    
    # approximate the rate integral and the survival probability integral for fast computations
    approxFuncs <- tess.prepare.pdf(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,max,t_crit)

    # compute the cumulative distribution function for the time of the process
    pdf <- function(x) tess.equations.pN.fastApprox(approxFuncs$r,approxFuncs$s,samplingProbability,nTaxa,0,x,SURVIVAL,MRCA,log=FALSE)

    
    obj.pdf <- function(t, state, pars) list(pdf(t))
    ## This step is slow for large n - perhaps 1/2s for 1000 points
    n1 <- 11
    times <- seq(0, max, length=n1)
    zz <- lsoda(0, times, obj.pdf, tcrit=max)[,2]
    m <- min(which( (zz[n1] - zz)/zz[n1] < 1E-5))
#    m2 <- max(which( zz/zz[n1] < 1E-5 )) # currently not used (Sebastian)
    if ( m > (n1/2) ) {
      # now we do it properly
      n2 <- 101
      times <- seq(0, max, length=n2)
      zz <- lsoda(0, times, obj.pdf, tcrit=max)[,2]
      zz <- zz / zz[n2]  # normalize
      break
    } else {
        max <- max/2.0
    }
  }

  icdf <- approxfun(zz, times) ## Interpolate
  simT <- function(n)  icdf(runif(n))
  T <- simT(n)
  
  trees <- list()
  # for each simulation
  for ( i in 1:n ) {

    # delegate the call to the simulation condition on both, age and nTaxa
    tree <- tess.sim.taxa.age.function(1,nTaxa,age=T[i],lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,samplingStrategy,MRCA,approxFuncs$r,approxFuncs$s)

    # add the new tree to the list
    trees[[i]] <- tree[[1]]
  }
  
  return (trees)
}

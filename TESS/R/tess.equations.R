################################################################################
#
# tess.equations.R
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
# @brief  Calculate the probability of survival in the interval [t_low,t_high].
#
# See Equation (2) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies under global, time-dependent birth-death processes
#
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t_low                                         scalar        starting time
# @param    t_high                                        scalar        end time
# @param    T                                             scalar        present time (time goes forward and the origin/MRCA might be 0)
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        probability of survival in [t,tau]
#
################################################################################
tess.equations.pSurvival.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log=FALSE) {

  # compute the rate
  rate <- mu - lambda

  # do the integration of int_{t_low}^{t_high} ( mu(s) exp(rate(t,s)) ds )
  # where rate(t,s) = int_{t}^{s} ( mu(x)-lambda(x) dx ) - sum_{for all t < m_i < s in massExtinctionTimes }( log(massExtinctionSurvivalProbability[i]) )

  # we compute the integral stepwise for each epoch between mass-extinction events
  # add mass-extinction
  accumulatedMassExtinction <- 1.0
  prev_time <- t_low
  den <- 1.0
  if ( length(massExtinctionTimes) > 0 ) {
    for (j in 1:length(massExtinctionTimes) ) {
      cond <-  (t_low < massExtinctionTimes[j]) & (t_high >= massExtinctionTimes[j])
      # compute the integral for this time episode until the mass-extinction event
#       den <- den + ifelse(cond, exp(-rate*t_low) * mu / (rate * accumulatedMassExtinction ) * ( exp(rate* massExtinctionTimes[j]) - exp(rate*prev_time)) , 0 )
      den <- den + cond * exp(-rate*t_low) * mu / (rate * accumulatedMassExtinction ) * ( exp(rate* massExtinctionTimes[j]) - exp(rate*prev_time))
      # store the current time so that we remember from which episode we need to integrate next
#       prev_time <- ifelse(cond, massExtinctionTimes[j], prev_time)
      prev_time[cond] <- massExtinctionTimes[j]
      accumulatedMassExtinction <- accumulatedMassExtinction * ifelse(cond, massExtinctionSurvivalProbabilities[j], 1.0)
      # integrate over the tiny time interval of the mass-extinction event itself and add it to the integral
#       den <- den - ifelse(cond, (massExtinctionSurvivalProbabilities[j]-1) / accumulatedMassExtinction * exp( rate*(massExtinctionTimes[j] - t_low) ), 0.0 )
      den <- den - cond * (massExtinctionSurvivalProbabilities[j]-1) / accumulatedMassExtinction * exp( rate*(massExtinctionTimes[j] - t_low) )
    }
  }

  # add the integral of the final epoch until the present time
  den <- den + exp(-rate*t_low) * mu / (rate * accumulatedMassExtinction ) * ( exp(rate*t_high) - exp(rate*prev_time))

  # add sampling
  cond <- (t_low < T) & (t_high >= T)
  if(samplingProbability < 1 ){
    accumulatedMassExtinction <- accumulatedMassExtinction * ifelse(cond, samplingProbability, 1.0)
  }

#   den <- den - ifelse(cond, (samplingProbability-1)*exp( rate*(T-t_low) ) / accumulatedMassExtinction, 0.0)
  den <- den - cond * (samplingProbability-1)*exp( rate*(T-t_low) ) / accumulatedMassExtinction

  res <- 1.0 / den

  if ( log == TRUE ) {
    res <- log( res )
  }

  return (res)

}



################################################################################
#
# @brief  Calculate the probability of survival in the interval [t_low,t_high].
#
# See Equation (2) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies
# under global, time-dependent birth-death processes
#
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        times at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t_low                                         scalar        starting time
# @param    t_high                                        scalar        end time
# @param    T                                             scalar        present time (time goes forward and the origin/MRCA might be 0)
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        probability of survival in [t,tau]
#
################################################################################
tess.equations.pSurvival.rateshift <- function(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log=FALSE,Rcpp=TRUE) {

  # do the integration of int_{t_low}^{t_high} ( mu(s) exp(rate(t,s)) ds )
  # where rate(t,s) = int_{t}^{s} ( mu(x)-lambda(x) dx ) - sum_{for all t < m_i < s in massExtinctionTimes }( log(massExtinctionSurvivalProbability[i]) )

  # we compute the integral stepwise for each epoch between mass-extinction events
  # add mass-extinction

  if(Rcpp){

  	if(is.null(rateChangeTimes)){
  		rateChangeTimes <- 0
  		massExtinctionSurvivalProbabilities <- 1
  	}

  	res <- equations_pSurvival_rateshift_CPP(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,t_low,t_high,T,log)
  	return(res)

  }

  accummulatedRateTime <- 0.0
  prev_time <- t_low
  den <- 1.0
  for (j in seq_len( length(rateChangeTimes) ) ) {
    # compute the rate
    rate <- mu[j] - lambda[j]
    cond <- (t_low < rateChangeTimes[j]) & (t_high >= rateChangeTimes[j])

    # compute the integral for this time episode until the mass-extinction event

    den <- den + cond*(exp(-rate*prev_time) * mu[j] / rate * exp( accummulatedRateTime ) * ( exp(rate* rateChangeTimes[j]) - exp(rate*prev_time)))
    accummulatedRateTime <- accummulatedRateTime + cond*(rate*(rateChangeTimes[j]-prev_time) - log(massExtinctionSurvivalProbabilities[j]))
    # store the current time so that we remember from which episode we need to integrate next
    prev_time[cond] <- rateChangeTimes[j]
    # integrate over the tiny time interval of the mass-extinction event itself and add it to the integral
    den <- den - cond*((massExtinctionSurvivalProbabilities[j]-1) * exp( accummulatedRateTime ))

  }

  if ( length(rateChangeTimes) > 0 ) {
    index <- min( length(lambda), findInterval(t_high,rateChangeTimes)+1)
  } else {
    index <- 1
  }
  rate <- mu[index] - lambda[index]

  # add the integral of the final epoch until the present time
  den <- den + exp(-rate*prev_time) * exp( accummulatedRateTime ) * mu[index] / rate * ( exp(rate*t_high) - exp(rate*prev_time))

  # add sampling
  cond <- (t_low < T) & (t_high >= T)
  accummulatedRateTime <- accummulatedRateTime + rate*(t_high-prev_time) - log(samplingProbability)
  den <- den - cond*((samplingProbability-1) * exp( accummulatedRateTime ))

  res <- 1.0 / den

  if ( log == TRUE ) {
    res <- log( res )
  }

  return (res)

}

################################################################################
#
# @brief  Calculate the probability of survival in the interval [t_low,t_high]
#         using a faster approximation.
#
# See Equation (2) from Hoehna, S., 2013, Fast simulation of reconstructed phylogenies
# under global, time-dependent birth-death processes
#
# @date Last modified: 2013-03-01
# @author Sebastian Hoehna
# @version 1.2
# @since 2012-09-18, version 1.0
#
# @param    r             function      rate integral function
# @param    s             function      survival probability function
# @param    t_low         scalar        starting time
# @param    t_high        scalar        end time
# @param    T             scalar        present time (time goes forward and the origin/MRCA might be 0)
# @param    log           bool          if in log-scale
# @return                 scalar        probability of survival in [t,tau]
#
################################################################################
tess.equations.pSurvival.fastApprox <- function(r, s, rho, t_low, t_high, T, log=FALSE) {

  if ( length(t_low) < length(t_high) )
    t_low <- rep( t_low, max(length(t_high),length(T)) )

  # Remember the following properties of r and s (assuming T = present):
  # r(x) = int_{0}{x} mu(t) - lambda(t) dt
  # s(x) = int_{x}{T} mu(t) exp( r(t) ) dt

  # The probability we want to compute is:
  # P( N(t_high)>0 | N(t_low)=1 ) = ( 1 + int_{t_low}^{t_high} mu(t)*exp(int_{t_low}^{t}mu(x)-lambda(x)dx) dt )^{-1}
  # and thus
  # P( N(t_high)>0 | N(t_low)=1 ) = ( 1 + exp(-r(t_low)) * (s(t_high)-s(t_low)) )^{-1}

  den <- ( 1 + ( ( s(t_low) - s(t_high) ) / exp(r(t_low)) ) )

  # add sampling
  cond <- (t_low < T) & (t_high >= T)
  den[cond] <- den[cond] - (rho-1) / rho *exp( r(T) - r(t_low) )[cond]
#  den[cond] <- den[cond] - (rho-1) / rho *exp( r(T) - r(t_low[cond]) ) # old code

  res <- ifelse( is.finite(den) & den > 0, 1/den, 0 )

  if ( log == TRUE )
    res <- log( res )

  return (res)
}


################################################################################
#
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 (or 2) lineage at time s.
#
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    i                                             scalar        number of lineages
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        log-probability
#
################################################################################
tess.equations.pN.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,i,s,t,SURVIVAL=FALSE,MRCA=FALSE,log=FALSE) {

  if ( i < 1 ) { # we assume conditioning on survival
    p <- 0
  } else if (i == 1) {
    if ( MRCA == TRUE ) { # we assume conditioning on survival of the two species
      p <- 0
    } else {
      if ( SURVIVAL == TRUE ) {
        p   <-  tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=TRUE) + (mu-lambda)*(t-s) - log(samplingProbability)
      } else {
        p   <-  2*tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=TRUE) + (mu-lambda)*(t-s) - log(samplingProbability)
      }
    }
  } else {
    p_s <- tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
    r   <- (mu-lambda)*(t-s) - log(samplingProbability)
    for (j in seq_len(length(massExtinctionTimes)) ) {
      cond <-  (s < massExtinctionTimes[j]) & (t >= massExtinctionTimes[j])
      r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
    }
    e   <- p_s * exp(r)
    e[e > 1] <- 1

    if ( MRCA == FALSE ) {
      if ( SURVIVAL == TRUE ) {
        p   <- log(p_s) + r + log( 1 - e) * (i-1)
      } else {
        p   <- 2*log(p_s) + r + log( 1 - e) * (i-1)
      }
    } else {
      if ( SURVIVAL == TRUE ) {
        p   <- log(i-1) + 2*log(p_s) + 2*r + log( 1 - e) * (i-2)
      } else {
        p   <- log(i-1) + 4*log(p_s) + 2*r + log( 1 - e) * (i-2)
      }
    }
  }

   if ( log == FALSE ) {
     p <- exp( p )
   }

  return (p)
}


################################################################################
#
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 (or 2) lineage at time s.
#
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    i                                             scalar        number of lineages
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        log-probability
#
################################################################################
tess.equations.pN.rateshift <- function(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,i,s,t,SURVIVAL=FALSE,MRCA=FALSE,log=FALSE) {

  # compute the rate
  prev_time <- t
  rate <- 0
  # add mass-extinction
  for (j in seq_len(length(rateChangeTimes)) ) {
    rate <- rate + ifelse( t < rateChangeTimes[j] & T >= rateChangeTimes[j], (mu[j] - lambda[j])*(rateChangeTimes[j]-prev_time) - log(massExtinctionSurvivalProbabilities[j]), 0 )
    prev_time <- ifelse( t < rateChangeTimes[j] & T >= rateChangeTimes[j], rateChangeTimes[j], t)
  }
  # add the final rate interval
  rate <- rate + ifelse( T > prev_time, (mu[length(mu)] - lambda[length(lambda)])*(T-prev_time), 0 )

  # add sampling
  rate <- rate - log(samplingProbability)


  if ( i < 1 ) { # we assume conditioning on survival
    p <- 0
  } else if (i == 1) {
    if ( MRCA == TRUE ) { # we assume conditioning on survival of the two species
      p <- 0
    } else {
      if ( SURVIVAL == TRUE ) {
        p   <-  tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=TRUE) + rate
      } else {
        p   <-  2*tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=TRUE) + rate
      }
    }
  } else {
    p_s <- tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
    e   <- p_s * exp(rate)
    e[e > 1] <- 1

    if ( MRCA == FALSE ) {
      if ( SURVIVAL == TRUE ) {
        p   <- log(p_s) + rate + log( 1 - e) * (i-1)
      } else {
        p   <- 2*log(p_s) + rate + log( 1 - e) * (i-1)
      }
    } else {
      if ( SURVIVAL == TRUE ) {
        p   <- log(i-1) + 2*log(p_s) + 2*rate + log( 1 - e) * (i-2)
      } else {
        p   <- log(i-1) + 4*log(p_s) + 2*rate + log( 1 - e) * (i-2)
      }
    }
  }

   if ( log == FALSE ) {
     p <- exp( p )
   }

  return (p)
}


################################################################################
#
# @brief  Calculate the probability of n lineage existing at time t
#         if we started with 1 (or 2) lineage(s) at time s. This is the fast
#         approximation.
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    rate                                          function      rate integral function
# @param    surv                                          function      survival integral function
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    i                                             scalar        number of lineages
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        log-probability
#
################################################################################
tess.equations.pN.fastApprox <- function(rate,surv,samplingProbability,i,s,t,SURVIVAL=FALSE,MRCA=FALSE,log=FALSE) {

  if ( i < 1 ) { # we assume conditioning on survival
    p <- 0
  } else if (i == 1) {
    if ( MRCA == TRUE ) { # we assume conditioning on survival of the two species
      p <- 0
    } else {
      if ( SURVIVAL == TRUE ) {
        p   <- tess.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=TRUE) + rate(t) - rate(s) - log(samplingProbability)
      } else {
        p   <- 2*tess.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=TRUE) + rate(t) - rate(s) - log(samplingProbability)
      }
    }
  }
  else {
    p_s <- tess.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=FALSE)
    r   <- rate(t) - rate(s) - log(samplingProbability)
    e   <- p_s * exp(r)
    e[e > 1] <- 1

    if ( MRCA == FALSE ) {
      if ( SURVIVAL == TRUE ) {
        p   <- log(p_s) + r + log( 1 - e) * (i-1)
      } else {
        p   <- 2*log(p_s) + r + log( 1 - e) * (i-1)
      }
    } else {
      if ( SURVIVAL == TRUE ) {
        p   <- log(i-1) + 2*log(p_s) + 2*r + log( 1 - e) * (i-2)
      } else {
        p   <- log(i-1) + 4*log(p_s) + 2*r + log( 1 - e) * (i-2)
      }
    }

  }

  if ( log == FALSE ) {
    p <- exp( p )
    p[is.finite(p) == FALSE] <- 0
  }


  return (p)
}


################################################################################
#
# @brief  Calculate the expected number of taxa when the process start at time
#         s with 1 (or two) species and ends at time t.
#
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @return                                                 scalar        expected number of taxa
#
################################################################################
tess.equations.nTaxa.expected.constant <- function(s,t,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA=FALSE) {

  p_s <- tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
  r   <- (mu-lambda)*(t-s) - log(samplingProbability)
  for (j in seq_len(length(massExtinctionTimes)) ) {
    cond <-  (s < massExtinctionTimes[j]) & (t >= massExtinctionTimes[j])
    r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
  }
  e   <- p_s * exp(r)
  e[e > 1] <- 1

  if ( MRCA == FALSE ) {
    n  <- 1.0 / e
  } else {
    n   <- 2.0 / e
  }

  return (n)
}


################################################################################
#
# @brief  Calculate the expected number of taxa when the process start at time
#         s with 1 (or two) species and ends at time t.
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    rate                                          function      rate integral function
# @param    surv                                          function      survival integral function
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @return                                                 scalar        the expected number of species
#
################################################################################
tess.equations.nTaxa.expected.fastApprox <- function(s,t,rate,surv,samplingProbability,MRCA=FALSE) {


  p_s <- tess.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=FALSE)
  r   <- rate(t) - rate(s) - log(samplingProbability)
  e   <- p_s * exp(r)
  e[e > 1] <- 1

  if ( MRCA == FALSE ) {
    n   <- 1.0 / e
  } else {
    n   <- 2.0 / e
  }

  return (n)
}



################################################################################
#
# @brief  Calculate the expected number of lineages at time t that will be
#         present in the reconstructed tree (i.e., non extinct),
#         when the process starts at time s with 1 (or two) species and ends at time 'present'.
#
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2014-02-11
# @author Sebastian Hoehna
# @version 2.0
# @since 2014-02-11, version 2.0
#
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @return                                                 scalar        expected number of taxa
#
################################################################################
tess.equations.nTaxa.expected.reconstructed.constant <- function(s,t,present,lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,MRCA=FALSE) {

  p_present <- tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,present,present,log=FALSE)
  p_s <- tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,s,t,t,log=FALSE)
  r   <- (mu-lambda)*(t-s) - log(samplingProbability)
  for (j in seq_len(length(massExtinctionTimes)) ) {
    cond <-  (s < massExtinctionTimes[j]) & (t >= massExtinctionTimes[j])
    r  <- r - ifelse(cond, log(massExtinctionSurvivalProbabilities[j]), 0.0)
  }
  e   <- ( 1.0   -   ( 1.0 - p_s * exp(r) )  * ( p_present / p_s )   )
  e[e > 1] <- 1

  if ( MRCA == FALSE ) {
    n  <- 1.0 / e
  } else {
    n   <- 2.0 / e
  }

  return (n)
}


################################################################################
#
# @brief  Calculate the expected number of taxa when the process start at time
#         s with 1 (or two) species and ends at time t.
#
# @see    Equation (3), (5) and (12) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-18, version 1.0
#
# @param    rate                                          function      rate integral function
# @param    surv                                          function      survival integral function
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    s                                             scalar        start time
# @param    t                                             scalar        stop time
# @param    MRCA                                          boolean       does the tree start at the mrca?
# @return                                                 scalar        the expected number of species
#
################################################################################
tess.equations.nTaxa.expected.reconstructed.fastApprox <- function(s,t,present,rate,surv,samplingProbability,MRCA=FALSE) {


  p_present <- tess.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,present,present,log=FALSE)
  p_s <- tess.equations.pSurvival.fastApprox(rate,surv,samplingProbability,s,t,t,log=FALSE)
  r   <- rate(t) - rate(s) - log(samplingProbability)
  e   <- ( 1.0   -   ( 1.0 - p_s * exp(r) )  * ( p_present / p_s )   )
  e[e > 1] <- 1

  if ( MRCA == FALSE ) {
    n   <- 1.0 / e
  } else {
    n   <- 2.0 / e
  }

  return (n)
}


################################################################################
#
# @brief  Calculate the probability of exactly 1 lineage surviving until time T
#         if we started with 1 lineage at time t.
#
# @see    Equation (4) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t                                             scalar        time
# @param    T                                             scalar        present time
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        ln-probability
#
################################################################################
tess.equations.p1.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=FALSE) {

  # compute the survival probability
  a <- tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=TRUE)

  # compute the rate
  rate <- (mu - lambda)*(T-t)
  # add mass-extinction
  for (j in seq_len(length(massExtinctionTimes)) ) {
    rate <- rate - ifelse( t < massExtinctionTimes[j] & T >= massExtinctionTimes[j], log(massExtinctionSurvivalProbabilities[j]), 0 )
  }

  # add sampling
  rate <- rate - log(samplingProbability)

  p <- 2*a + rate

  if ( log == FALSE ) {
    p <- exp( p )
  }

  return (p)
}


################################################################################
#
# @brief  Calculate the probability of exactly 1 lineage surviving until time T
#         if we started with 1 lineage at time t.
#
# @see    Equation (4) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    lambda                                        scalar        speciation rate
# @param    mu                                            scalar        extinction rate
# @param    massExtinctionTimes                           vector        timse at which mass-extinctions happen
# @param    massExtinctionSurvivalProbabilities           vector        survival probability of a mass extinction event
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t                                             scalar        time
# @param    T                                             scalar        present time
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        ln-probability
#
################################################################################
tess.equations.p1.rateshift <- function(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=FALSE,Rcpp=TRUE) {

  if(Rcpp){

  	if(is.null(rateChangeTimes)){
  		rateChangeTimes <- 0
  		massExtinctionSurvivalProbabilities <- 1
  	}

  	p <- equations_p1_rateshift_CPP(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,log=TRUE)

  	return(p)

  }

  # compute the survival probability
  a <- tess.equations.pSurvival.rateshift(lambda,mu,rateChangeTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=TRUE,Rcpp=Rcpp)

  # compute the rate
  prev_time <- t
  rate <- 0
  # add mass-extinction
  for (j in seq_len(length(rateChangeTimes)) ) {

 	true <- t < rateChangeTimes[j] & T >= rateChangeTimes[j]
 	rate <- rate + true*((mu[j] - lambda[j])*(rateChangeTimes[j]-prev_time) - log(massExtinctionSurvivalProbabilities[j]))
 	prev_time[true] <- rateChangeTimes[j]

  }

  rate <- rate + (T > prev_time)*((mu[length(mu)] - lambda[length(lambda)])*(T-prev_time))

  # add sampling
  rate <- rate - log(samplingProbability)

  p <- 2*a + rate

  if ( log == FALSE ) {
    p <- exp( p )
  }

  return (p)
}


################################################################################
#
# @brief  Calculate the probability of exactly 1 lineage surviving until time T
#         if we started with 1 lineage at time t.
#
# @see    Equation (4) in Hoehna, S.: Fast simulation of reconstructed phylogenies
#         under global, time-dependent birth-death processes. 2013, Bioinformatics
#
# @date Last modified: 2013-01-30
# @author Sebastian Hoehna
# @version 1.1
# @since 2012-09-11, version 1.0
#
# @param    r                                             function      rate integral function
# @param    s                                             function      survival integral function
# @param    samplingProbability                           scalar        probability of uniform sampling at present
# @param    t                                             scalar        time
# @param    T                                             scalar        present time
# @param    log                                           bool          if in log-scale
# @return                                                 scalar        ln-probability
#
################################################################################
tess.equations.p1.fastApprox <- function(r,s,samplingProbability,t,T,log=FALSE) {
  # Remember the following properties of r and s (assuming T = present):
  # r(x) = int_{0}{x} mu(t) - lambda(t) dt
  # s(x) = int_{x}{T} mu(t) exp( r(t) ) dt

  # The probability we want to compute is:
  # P( N(T)=1 | N(t)=1 ) = P( N(T)>0 | N(t)=1 )^{2} * exp( int_{t}^{T} mu(x)-lambda(x) dx )
  # and thus
  # P( N(t_high)>0 | N(t_low)=1 ) = ( 1 + exp(-r(t_low)) * (s(t_low)-s(t_high)) )^{-1}

  # for simplicity we compute the log-probability
  a <- tess.equations.pSurvival.fastApprox(r,s,samplingProbability,t,T,T,log=TRUE)
  b <- r(T) - r(t) - log(samplingProbability) # Note that we always include the present time and thus always apply uniform taxon sampling
  p <- 2*a + b

  if ( log == FALSE ) {
    p <- exp( p )
  }

  return (p)
}


################################################################################
#
# @brief  Calculate the probability of no speciation event on the reconstructed
#         process between time t and t'.
#
# @date Last modified: 2013-02-06
# @author Sebastian Hoehna
# @version 1.2
# @since 2013-02-06, version 1.2
#
# @param    lambda        scalar        speciation rate
# @param    mu            scalar        extinction rate
# @param    t             scalar        time
# @param    T             scalar        present time
# @param    log           bool          if in log-scale
# @return                 scalar        ln-probability
#
################################################################################
tess.equations.pWaiting.constant <- function(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,t_prime,T,n,log=FALSE) {

  u <- 1.0 - tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,t_prime,T,log=FALSE) * exp( (mu-lambda)*(t_prime-t) - log(samplingProbability) )
  tmp <- u * tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,T,T,log=FALSE) / tess.equations.pSurvival.constant(lambda,mu,massExtinctionTimes,massExtinctionSurvivalProbabilities,samplingProbability,t,t_prime,T,log=FALSE)

  # some modification so that we can compute the log
  tmp[tmp > 1.0] <- 1.0
  p <- log(1.0 - tmp ) * n

  if (log == FALSE) {
    p <- exp(p)
  }

  return(p)

}





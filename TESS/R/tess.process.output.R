################################################################################
#
# tess.plot.output.R
#
# Copyright (c) 2012- Michael R May
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
# @brief Processing the output of a episodic diversification rate analysis with mass-extinction events.
#
# @date Last modified: 2014-10-05
# @author Michael R May
# @version 2.0
# @since 2014-10-04, version 2.0.0
#
# @param    dir                         character      The directory from which the CoMET output will be read.
# @param    tree                        phylo          The tree analyzed with CoMET in phylo format. By default, looks for a tree in the target directory.
# @param    numExpectedRateChanges      numeric        The number of expected diversification-rate changes.
# @param    numExpectedMassExtinctions  numeric        The number of expected mass-extinction events.
# @param    burnin                      numeric        The fraction of the posterior samples to be discarded as burnin.
# @param    numIntervals                numeric        The number of discrete intervals in which to break the tree.
# @param    criticalBayesFactors        numeric        The Bayes factor thresholds to use to assess significance of events.
#
#
################################################################################

tess.process.output = function(dir,tree=NULL,numExpectedRateChanges=2,numExpectedMassExtinctions=2,burnin=0.25,numIntervals=100,criticalBayesFactors=c(2,6,10)){

  files <- list.files(dir,full.names=TRUE)

  if ( !any(grepl("samples_numCategories",files)) ) {
    stop("The directory provided does not appear to contain any CoMET analysis output.")
  }

  if( is.null(tree) ) {
    if ( any(grepl(".tre",list.files(dir))) ) {
      tree <- read.nexus(grep(".tre",files,fixed=TRUE,value=TRUE))
    } else {
      stop("You must either provide an input tree or have a correct .tre file into the directory.")
    }
  }

  # Get the time of the tree and divide it into intervals
  time <- max( branching.times(tree) )
  intervals <- seq(0,time,length.out=numIntervals+1)

  # Get the numCategories output
  numCategoriesOutput <- read.table(grep("samples_numCategories",files,value=TRUE),header=TRUE)
  categoriesBurnin <- nrow(numCategoriesOutput) * burnin

  # Process the posterior probality of the model
  modelPosteriorProbability <- as.mcmc(numCategoriesOutput$posterior[categoriesBurnin:nrow(numCategoriesOutput)])

  # Process the number of speciation rate categories
  speciationRateCategories <- as.mcmc(numCategoriesOutput$NumSpeciation[categoriesBurnin:nrow(numCategoriesOutput)])

  # Process the number of extinction rate categories
  extinctionRateCategories <- as.mcmc(numCategoriesOutput$numExtinction[categoriesBurnin:nrow(numCategoriesOutput)])

  # Process the number of mass-extinction events
  numMassExtinctions <- as.mcmc(numCategoriesOutput$numMassExtinctions[categoriesBurnin:nrow(numCategoriesOutput)])

  # Process the speciation rates
  speciationRateChangeTimes <- strsplit(readLines(grep("SpeciationRateChanges",files,value=TRUE))[-1],"\t")
  speciationRates <- strsplit(readLines(grep("SpeciationRates",files,value=TRUE))[-1],"\t")
  speciationBurnin <- length(speciationRates) * burnin

  processSpeciationRates <- as.mcmc(do.call(rbind,lapply(speciationBurnin:length(speciationRateChangeTimes),function(sample) {
    times <- as.numeric(speciationRateChangeTimes[[sample]])
    rates <- as.numeric(speciationRates[[sample]])
    order <- order(times)
    times <- times[order]
    rates <- c(rates[1],rates[-1][order])
    res   <- rates[findInterval(intervals[-1],times)+1]
    return (res)
  } )))

  processSpeciationRateChangeTimes <- as.mcmc(do.call(rbind,lapply(1:length(speciationRateChangeTimes),function(sample) {
    times <- sort(as.numeric(speciationRateChangeTimes[[sample]]))
    res <- rep(0,numIntervals)
    res[findInterval(times,intervals)] <- 1
    return (res)
  } )))

  # Compute the Bayes factors for speciation rate change events
  speciationRateChangePosteriorProbability <- colMeans(processSpeciationRateChangeTimes)
  speciationRateChangePriorProbability     <- 1 - dpois(0,lambda=numExpectedRateChanges/numIntervals)
  speciationRateChangePosteriorModelOdds   <- ( speciationRateChangePosteriorProbability / (1 - speciationRateChangePosteriorProbability) )
  speciationRateChangePriorModelOdds       <- ( speciationRateChangePriorProbability / (1 - speciationRateChangePriorProbability) )
  speciationRateChangeBayesFactors         <- 2 * log( speciationRateChangePosteriorModelOdds / speciationRateChangePriorModelOdds )

  # Process the extinction rates
  extinctionRateChangeTimes <- strsplit(readLines(grep("ExtinctionRateChanges",files,value=TRUE))[-1],"\t")
  extinctionRates <- strsplit(readLines(grep("ExtinctionRates",files,value=TRUE))[-1],"\t")
  extinctionBurnin <- length(extinctionRates) * burnin

  processExtinctionRates <- as.mcmc(do.call(rbind,lapply(extinctionBurnin:length(extinctionRateChangeTimes),function(sample) {
    times <- as.numeric(extinctionRateChangeTimes[[sample]])
    rates <- as.numeric(extinctionRates[[sample]])
    order <- order(times)
    times <- times[order]
    rates <- c(rates[1],rates[-1][order])
    res   <- rates[findInterval(intervals[-1],times)+1]
    return (res)
  } )))

  processExtinctionRateChangeTimes <- as.mcmc(do.call(rbind,lapply(1:length(extinctionRateChangeTimes),function(sample) {
    times <- sort(as.numeric(extinctionRateChangeTimes[[sample]]))
    res <- rep(0,numIntervals)
    res[findInterval(times,intervals)] <- 1
    return (res)
  } )))

  # Compute the Bayes factors for extinction rate change events
  extinctionRateChangePosteriorProbability <- colMeans(processExtinctionRateChangeTimes)
  extinctionRateChangePriorProbability     <- 1 - dpois(0,lambda=numExpectedRateChanges/numIntervals)
  extinctionRateChangePosteriorModelOdds   <- ( extinctionRateChangePosteriorProbability / (1 - extinctionRateChangePosteriorProbability) )
  extinctionRateChangePriorModelOdds       <- ( extinctionRateChangePriorProbability / (1 - extinctionRateChangePriorProbability) )
  extinctionRateChangeBayesFactors         <- 2 * log( extinctionRateChangePosteriorModelOdds / extinctionRateChangePriorModelOdds )

  # Process the net-diversification and relative-extinction rates
  processNetDiversificationRates <- as.mcmc(processSpeciationRates-processExtinctionRates)
  processRelativeExtinctionRates <- as.mcmc(processExtinctionRates/processSpeciationRates)

  # Process the mass extinctions
  massExtinctionTimes <- strsplit(readLines(grep("MassExtinctionTimes",files,value=TRUE))[-1],"\t")
  massExtinctionBurnin <- length(massExtinctionTimes) * burnin

  processMassExtinctionTimes <- as.mcmc(do.call(rbind,lapply(massExtinctionBurnin:length(massExtinctionTimes),function(sample) {
    times <- sort(as.numeric(massExtinctionTimes[[sample]]))
    res <- rep(0,numIntervals)
    res[findInterval(times,intervals)] <- 1
    return (res)
  } )))

  # Compute the Bayes factors for mass extinction events
  massExtinctionPosteriorProbability <- colMeans(processMassExtinctionTimes)
  massExtinctionPriorProbability     <- 1 - dpois(0,lambda=numExpectedMassExtinctions/numIntervals)
  massExtinctionPosteriorModelOdds   <- ( massExtinctionPosteriorProbability / (1 - massExtinctionPosteriorProbability) )
  massExtinctionPriorModelOdds       <- ( massExtinctionPriorProbability / (1 - massExtinctionPriorProbability) )
  massExtinctionBayesFactors         <- 2 * log( massExtinctionPosteriorModelOdds / massExtinctionPriorModelOdds )

  # Compute the critical Bayes factors
  criticalBayesFactors <- exp( c(2,6,10) / 2 )
  x <- criticalBayesFactors * massExtinctionPriorModelOdds
  massExtinctionCriticalPosteriorProbabilities <- x / (1 + x)
  x <- criticalBayesFactors * extinctionRateChangePriorModelOdds
  speciationRateChangeCriticalPosteriorProbabilities <-
  extinctionRateChangeCriticalPosteriorProbabilities <- x / (1 + x)

  res <- list("posterior" = modelPosteriorProbability,
              "numSpeciationCategories" = speciationRateCategories,
              "numExtinctionCategories" = extinctionRateCategories,
              "numMassExtinctions" = numMassExtinctions,
              "speciation rates" = processSpeciationRates,
              "speciation shift times" = processSpeciationRateChangeTimes,
              "speciation Bayes factors" = speciationRateChangeBayesFactors,
              "speciationRateChangeCriticalPosteriorProbabilities" = speciationRateChangeCriticalPosteriorProbabilities,
              "extinction rates" = processExtinctionRates,
              "extinction shift times" = processExtinctionRateChangeTimes,
              "extinction Bayes factors" = extinctionRateChangeBayesFactors,
              "extinctionRateChangeCriticalPosteriorProbabilities" = extinctionRateChangeCriticalPosteriorProbabilities,
              "net-diversification rates" = processNetDiversificationRates,
              "relative-extinction rates" = processRelativeExtinctionRates,
              "mass extinction times" = processMassExtinctionTimes,
              "mass extinction Bayes factors" = massExtinctionBayesFactors,
              "massExtinctionCriticalPosteriorProbabilities" = massExtinctionCriticalPosteriorProbabilities,
              "criticalBayesFactors" = criticalBayesFactors,
              "tree" = tree,
              "intervals" = rev(intervals) )

  return(res)

}


































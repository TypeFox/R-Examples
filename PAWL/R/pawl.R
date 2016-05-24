###################################################
#    This file is part of RPAWL.
#
#    RPAWL is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RPAWL is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RPAWL.  If not, see <http://www.gnu.org/licenses/>.
###################################################
## Function to check if the flat histogram criterion is reached
checkFlatHistogram <- function(FHbincount, binning){
  difference <- abs((FHbincount / sum(FHbincount)) - binning@desiredfreq)
  FHreached <- all(difference < binning@fhthreshold * binning@desiredfreq) &
  all(FHbincount > 10)
  return(FHreached)
}


## Particle Wang-Landau function with flat histogram criterion
## AP stands for Algorithmic Parameters
pawl <- function(target, binning, AP, proposal, verbose = TRUE){
    if (verbose) print("Launching Particle Wang-Landau algorithm ...") 
    # Init some algorithmic parameters ...
    nbins <- length(binning@bins)
    # The bias is saved in a matrix but because
    # the dimension might change (if bins are split)
    # we also store a list of matrices of possibly different
    # dimensions, if bins have been split.
    logtheta <- matrix(ncol = nbins, nrow = AP@niterations + 1)
    logtheta[1,] <- 0
    logthetahistory <- list()
    # keep track of the count of visits to each bin
    # since the last flat histogram
    if (binning@useFH)
        FHbincount <- rep(0, nbins)
    # We compute the bins' middles in case the automatic bin split
    # mechanism is enabled. We also keep the count of visits in the bins,
    # as well as the count of visits in the half left part of the bins ...
    binning@binmids <- getBinMiddles(binning)
    diagnoseactive <- binning@autobinning
    if (diagnoseactive){
        diagnosebincount <- rep(0, nbins)
        diagnosehalfleftcount <- rep(0, nbins - 2)
    }
    # k is the temperature of the stochastic approximation tempering,
    # that makes the update of the bias less and less important. k is
    # increased when the flat histogram criterion is met.
    k <- 1
    khistory <- c(k)
    FHtimes <- c()
    # Due to the bin split mechanism, we keep track of the bins ...
    splitTimes <- c()
    binshistory <- list()
    binshistory[[1]] <- binning@bins
    nbinsvector <- c(nbins)
    # We keep track of the chains, and of the history in "allchains"...
    chains <- as.matrix(target@rinit(AP@nchains))
    # We keep track of the log densities computed along the iterations ...
    currentlogtarget <- target@logdensity(chains, target@parameters)
    if (AP@computemean){
        sumchains <- matrix(0, nrow = AP@nchains, ncol = target@dimension)
    }
    acceptrates <- rep(0, AP@niterations)
    # Setting the proposal distribution for the MH kernel
    if (missing(proposal) & target@type == "continuous"){
        proposal <- createAdaptiveRandomWalkProposal(nchains = AP@nchains, 
                                               targetdimension = target@dimension,
                                               adaptiveproposal = TRUE)
    }
    if (proposal@adaptiveproposal){
        sigma <- rep(0, AP@niterations + 1)
        sigma[1] <- proposal@proposalparam$sigma
    }
    dproposal <- proposal@dproposal
    rproposal <- proposal@rproposal
    proposalparam <- proposal@proposalparam
    # standard dev of the adaptive proposal ...
    if (proposal@adaptiveproposal){
        sigma <- rep(0, AP@niterations + 1)
        sigma[1] <- proposal@proposalparam$sigma
    }
    if (target@type == "discrete" & proposal@adaptiveproposal){
        proposal@adaptiveproposal <- FALSE
        if (verbose) cat("switching off adaptive proposal, because the target is discrete\n")
    }
    # We compute the locations of the chains, that is, in which 
    # bins they are.
    currentreaction <- binning@position(chains, currentlogtarget)
    currentlocations <- binning@getLocations(binning@bins, currentreaction)
    thisiterationcount <- tabulate(currentlocations, nbins = nbins)
    if (binning@useFH)
        FHbincount <- FHbincount + thisiterationcount
    if (diagnoseactive){
      diagnosebincount <- thisiterationcount
      diagnosehalfleftcount <- getInnerLeftCounts(binning,
                                                  currentreaction, 
                                                  currentlocations)
    }
    lastFHtime <- 0
    lastSplittime <- 0
    if (AP@saveeverynth > 0){
      nallchains <- 1 + floor((AP@niterations) / AP@saveeverynth)
      if (verbose) cat("saving every", AP@saveeverynth, "iterations\n")
      totalsize <- nallchains * AP@nchains * target@dimension
      if (verbose)
          cat("hence saving a vector of size", nallchains, "x", AP@nchains, 
              "x", target@dimension, "=", totalsize, "\n")
      if (totalsize > 10^8){
        if (verbose) cat("which bigger than 10^8: you better have a lot of memory available!!\n")
        suggestedmaxnallchains <- floor((10^8) / (target@dimension * AP@nchains)) + 1
        if (verbose) cat("you can maybe set saveeverynth to something bigger than ", 
            floor(AP@niterations / suggestedmaxnallchains) + 1, "\n")
        if (verbose) cat("type anything to continue, or Ctrl-C to abort\n")
        y<-scan(n=1)
      }
      allchains <- array(NA, dim = c(nallchains, AP@nchains, target@dimension))
      nstoredchains <- 1
      allchains[nstoredchains,,] <- chains
      alllogtarget <- matrix(NA, nrow = nallchains, ncol = AP@nchains)
      alllogtarget[nstoredchains,] <- currentlogtarget
    }
    allreaction <- matrix(nrow = AP@niterations + 1, ncol = AP@nchains)
    allreaction[1,] <- currentreaction
    iterstep <- max(100, AP@niterations / 50)
    for (iteration in 1:AP@niterations){
        if (!(iteration %% iterstep)){
          if (verbose) cat("iteration", iteration, "/", AP@niterations, "\n")
        }
        ## Sample new values from the MH kernel targeting the biased distribution ...
        rproposalresults <- rproposal(chains, proposalparam)
        proposals <- rproposalresults$states
        if (target@updateavailable){
            proposalLogTarget <- currentlogtarget + target@logdensityupdate(chains, 
                                              target@parameters, rproposalresults$others)
        } else {
            proposalLogTarget <- target@logdensity(proposals, target@parameters)
        }
        proposalReaction <- binning@position(proposals, proposalLogTarget)
        proposalLocations <- binning@getLocations(binning@bins, proposalReaction)
        loguniforms <- log(runif(AP@nchains))

        accepts <- (loguniforms < ((proposalLogTarget + dproposal(proposals, chains, proposalparam)
                                    - logtheta[iteration - lastSplittime,][proposalLocations]) 
        - (currentlogtarget + dproposal(chains, proposals, proposalparam)
           - logtheta[iteration - lastSplittime,][currentlocations ])))
        chains[accepts,] <- proposals[accepts,]
        currentlogtarget[accepts] <- proposalLogTarget[accepts]
        currentlocations[accepts] <- proposalLocations[accepts]
        currentreaction[accepts] <- proposalReaction[accepts]
        if (AP@saveeverynth > 0 & iteration %% AP@saveeverynth == 0){
            nstoredchains <- nstoredchains + 1
            allchains[nstoredchains,,] <- chains
            alllogtarget[nstoredchains,] <- currentlogtarget
        }
        allreaction[iteration + 1,] <- currentreaction
        if (AP@computemean && iteration > AP@computemeanburnin){
            sumchains <- sumchains + chains
        }
        acceptrates[iteration] <- mean(accepts)
        ## update the proportions of visit in each bin
        ## and the inner distribution of the chains in each bin
        ## (proportions in the left hand side of each bin)
        thisiterationcount <- tabulate(currentlocations, nbins = nbins)
        if (binning@useFH)
            FHbincount <- FHbincount + thisiterationcount
        if (diagnoseactive){
            diagnosebincount <- diagnosebincount + thisiterationcount
            diagnosehalfleftcount <- diagnosehalfleftcount + getInnerLeftCounts(binning, currentreaction, currentlocations)
        }
        ## update the bias using all the chains ...
        currentproportions <- thisiterationcount / AP@nchains
        if (binning@useLearningRate){
            if (binning@useFH){
                logtheta[iteration + 1 - lastSplittime,] <- logtheta[iteration - 
                             lastSplittime,] + binning@learningrate(k) * 
                             (currentproportions - binning@desiredfreq)
            } else {
                logtheta[iteration + 1 - lastSplittime,] <- logtheta[iteration - 
                             lastSplittime,] + binning@learningrate(iteration) * 
                             (currentproportions - binning@desiredfreq)
            }
        } else {
            logtheta[iteration + 1 - lastSplittime,] <- logtheta[iteration - lastSplittime,] + 
                    (currentproportions - binning@desiredfreq)
        }
        ## update the adaptive proposal standard deviation ...
        if (proposal@adaptiveproposal){
            sigma[iteration + 1] <- max(10^(-10 - target@dimension), sigma[iteration] + 
                proposal@adaptationrate(iteration) * (2 * (mean(accepts) > 0.234) - 1))
            proposalparam$sigma <- sigma[iteration + 1]
        } 
        if (diagnoseactive & binning@alongenergy){
            if (diagnosebincount[nbins] > AP@nchains){
                if (verbose) cat("right end bin reached: disactivating bin diagnosis\n")
                diagnoseactive <- FALSE
            }
        }
        if (diagnoseactive & (iteration %% max(100, 1 + floor(AP@niterations / 100)) == 0) &
            iteration < AP@niterations){
            if (verbose) cat("/** diagnosis at iteration", iteration, "\n")
            if (verbose) cat("* current number of bins =", nbins, "\n")
            allproportions <- diagnosebincount / sum(diagnosebincount)
            foundbins <- findSkewedBins(binning, 
                                         diagnosehalfleftcount, diagnosebincount)
            skewedbins <- foundbins$skewedbins + 1
            problem_bins <- (binning@desiredfreq - allproportions) / binning@desiredfreq
            if (verbose) cat("desired freq - proportions:\n", (problem_bins), "\n")
            if (verbose) cat("* skewed bins:", skewedbins, "\n")
            if (verbose) cat("* bins with enough points to be split:", 
                which(diagnosebincount > 10 * AP@nchains), "\n")
            bintosplit <- c()
            if (binning@alongenergy){
                for (binindex in 2:(nbins - 1)){
                    if (diagnosebincount[binindex] > 10 * AP@nchains){
                        if (all(problem_bins[(binindex + 1):nbins] > 0.9)){
                            bintosplit <- c(bintosplit, binindex)
                        }
                    }  
                }
            } else {
                for (binindex in 2:(nbins - 1)){
                    if (diagnosebincount[binindex] > 10 * AP@nchains &
                        binindex %in% skewedbins){
                        if (binning@smoothbinning | problem_bins[binindex + 1] > 0.9 |
                            problem_bins[binindex - 1] > 0.9){
                            bintosplit <- c(bintosplit, binindex)
                        }
                    }  
                }
            }
            if (!is.null(bintosplit)){
                if (verbose) cat("* bins to split:", bintosplit, "\n")
                newcuts <- binning@binmids[bintosplit - 1]
                foundbins <- list(binsToSplit = bintosplit, newcuts = newcuts)
                splitresults <- binsplitter(binning, foundbins, 
                    logtheta[iteration + 1 - lastSplittime,], binning@desiredfreq, 
                    binning@alongenergy)
                newbins <- splitresults$newbins
                nbins <- length(newbins)
                binning@bins <- newbins
                binning@binmids <- getBinMiddles(binning)
                binning@desiredfreq <- splitresults$newdesiredfreq
                if (binning@useFH){
                    k <- 1
                    khistory[length(khistory)] <- 1
                    FHbincount <- rep(0, nbins)
                }
                diagnosebincount <- rep(0, nbins)
                diagnosehalfleftcount <- rep(0, nbins - 2)
                # update the current values that depend on bins
                currentreaction <- binning@position(chains, currentlogtarget)
                currentlocations <- binning@getLocations(binning@bins, 
                                                         currentreaction)
                # store things
                nbinsvector <- c(nbinsvector, nbins)
                binshistory[[length(binshistory) + 1]] <- newbins
                lastSplittime <- iteration
                splitTimes <- c(splitTimes, iteration)
                if (length(splitTimes) == 1){
                    logthetahistory[[1]] <- logtheta[1:(iteration + 1),]
                } else {
                    diffSplitTimes <- splitTimes[length(splitTimes)] - 
                        splitTimes[length(splitTimes) - 1]
                    logthetahistory[[length(splitTimes)]] <- 
                        logtheta[1:(diffSplitTimes + 1),]
                }
                logtheta <- matrix(ncol = nbins, nrow = AP@niterations + 1 - iteration)
                logtheta[1,] <- splitresults$newthetas
            }
            if (verbose) cat("*/\n")
        }

        ## check if flat histogram is reached ...
        if (binning@useFH){
            FHreached <- checkFlatHistogram(FHbincount, binning)
            enoughElapsedTime <- (iteration >= lastFHtime + binning@minSimEffort)
            if (FHreached && enoughElapsedTime){
                if (verbose) cat("Flat histogram criterion met at iteration", iteration, "!\n")
                k <- k + 1
                khistory <- c(khistory, k)
                FHtimes <- c(FHtimes, iteration)
                lastFHtime <- iteration
                if (!binning@smoothbinning)
                  diagnoseactive <- FALSE
                FHbincount <- rep(0, nbins)
            }
        }
    }
    results <- list(chains = chains, acceptrates = acceptrates, logtheta = logtheta,
                finallocations = currentlocations, FHtimes = FHtimes, 
                finalbins = binning@bins, finaldesiredfreq = binning@desiredfreq,
                splitTimes = splitTimes, nbins = nbinsvector,
                binshistory = binshistory, khistory = khistory)
    if (proposal@adaptiveproposal)
        results$sigma <- sigma
    if (AP@saveeverynth > 0){
        results$allchains <- allchains
        results$alllogtarget <- alllogtarget
    }
    results$allreaction <- allreaction
    if (AP@computemean){
        results$meanchains <- sumchains / (AP@niterations - AP@computemeanburnin)
    }
    logthetahistory[[length(logthetahistory) + 1]] <- logtheta
    results$logthetahistory <- logthetahistory
    return(results)
}

getFrequencies <- function(results, binning){
    allproportions <- tabulate(binning@getLocations(results$finalbins, 
                      c(results$allreaction)), nbins = length(results$finalbins))
    finalfrequencies <- allproportions / sum(allproportions)
    innerfinalbins <- results$finalbins
    innerfinalbins <- innerfinalbins[2:(length(innerfinalbins))]
    innerinitbins <- results$binshistory[[1]]
    innerinitbins <- innerinitbins[2:(length(innerinitbins))]
    samplefrequencies <- c(finalfrequencies[1])
    for (index in 1:(length(innerinitbins)-1)){
        indexstart <- which(innerinitbins[index] == innerfinalbins)
        indexstop <- which(innerinitbins[index+1] == innerfinalbins)
        samplefrequencies <- c(samplefrequencies, sum(finalfrequencies[(1+indexstart):indexstop]))
    }
    samplefrequencies <- c(samplefrequencies, finalfrequencies[length(finalfrequencies)])
    cat("/* Frequency check\n")
    cat("*Do the obtained frequencies match the desired frequencies?\n")
    cat("*final bins:", results$finalbins, "\n")
    cat("*corresponding frequencies:", finalfrequencies, "\n")
    cat("*initial bins:", results$binshistory[[1]], "\n")
    cat("*desired frequencies: ", binning@desiredfreq, "\n")
    cat("*obtained frequencies:", samplefrequencies, "\n")
    return(samplefrequencies)
}




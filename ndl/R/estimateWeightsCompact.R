#### Companion function to estimateWeights, but uses compact binary event data files
#### See R documentation for details.

estimateWeightsCompact <- 
  function(datasource, removeDuplicates=TRUE, saveCounts=FALSE, verbose=FALSE, MaxEvents=100000000000000, trueCondProb=TRUE, addBackground=FALSE, ...)
{
    
    ## Check for pre-existing counts.
    coocFile = paste(datasource,".coocCues.rds",sep='')
    coocOutFile = paste(datasource,".coocCuesOutcomes.rds",sep='')
    if ((file.exists(coocFile) & file.exists(coocOutFile)) & MaxEvents >= 100000000000000 ) {
        warning("Did not remove duplicates because there were pre-computed cooccurrence matrices availbe. Remove these files and run again.")
        message(paste(c("\nNOTE: Loading pre-computed coocurrence matrices.\nIgnoring DataFrame '", basename, "' Provided.\nPlease remove the files ",coocFile," and ",coocOutFile, " if this behavior is not desired.")),sep="")
        flush.console()
        coocCues= readRDS(coocFile)
        coocCuesOutcomes= readRDS(coocOutFile)
    } else {
        if (verbose) message("Reading compact binary data from disk.")
        flush.console()
        ## call external C++ function using RCpp
        CuAndCo = learn(data=datasource,RemoveDuplicates=removeDuplicates, verbose=verbose, MaxEvents=MaxEvents, addBackground=addBackground) 
        coocCues = CuAndCo[[1]]
        coocCuesOutcomes = CuAndCo[[2]]
        coocOutcomesFreq = CuAndCo[[3]]
        if ((nrow(coocCuesOutcomes) <2) | (ncol(coocCuesOutcomes) <2)) {
            stop("Your data had insufficient number of unique cues or outcomes. Please make sure that you have at least two cues and at least two outcomes.")
        }
        ## Save the cooc matrices for later reuse (after doing Background rates and normalization.
        if (saveCounts) {
            if (verbose) message("Completed Event Counts. Saving so-occurrence data for future calculations.")
            flush.console()
            saveRDS(coocCues, file=coocFile)
            if (verbose) message(paste("Saved",coocFile))
            flush.console()
            saveRDS(coocCuesOutcomes, file=coocOutFile)
            if (verbose) message(paste("Saved",coocOutFile))
            flush.console()
        }
    }
    if (verbose) message("Starting to process matrices.")
    ## Check sanity of arguments
    if ((addBackground) & (!trueCondProb)) {
        message("*WARNING: Can't add background rates without true conditional probabilities. \n*ACTION: Proceeding without background rates.")
        addBackground = FALSE
    }
    if (addBackground & trueCondProb) {
      # Add background for Cue-Cue cooc
      Environ = diag(coocCues)
      grandTotal = sum(Environ)
      coocCues = rbind(Environ,coocCues)
      Environ=c(grandTotal,Environ)
      coocCues = cbind(Environ,coocCues)
      rownames(coocCues)[1]=c("Environ")
      # Add background for Cue-Outcome cooc
      Environ=coocOutcomesFreq 
      coocCuesOutcomes = rbind(Environ,coocCuesOutcomes)
    }
    if (trueCondProb) {
      ## Convert Cue-Outcome counts to Cue-Outcome Probabilities using diagonal
      cueTotals = diag(coocCues) 
      cueTotals[cueTotals == 0] = 1
      condProbsCues = coocCues/cueTotals
      probsOutcomesGivenCues = coocCuesOutcomes/cueTotals
    } else {
      ## use the original algorithm for normalization
      rowsums = rowSums(coocCuesOutcomes)
      rowsums[rowsums == 0] = 1
      condProbsCues = coocCues/rowsums
      probsOutcomesGivenCues = coocCuesOutcomes/rowsums
    }

    if (verbose) message("Starting to calculate pseudoinverse.")
    flush.console()
    n = dim(condProbsCues)[1]
    ## Check to see if the number of cues is too big for reasonable hardware.
    if (n < 20000) {
      pseudoinverse = ginv(condProbsCues)
    } else {
        ## Use an approximation of the pseudoinverse here to make this feasible
        ## average hardware.
        if (verbose) message("Number of cues was too large for standard pseudoinverse. Switching to lower-rank approximation.")
        pseudoinverse = random.pseudoinverse(condProbsCues,verbose=verbose)
    }
    ## Calculate the weights by multiplying the pseudoinver of the c-c
    ## counts by the probabilites of the outcomes given the cues.
    weightMatrix = pseudoinverse %*% probsOutcomesGivenCues
    rownames(weightMatrix) = rownames(coocCues)
    colnames(weightMatrix) = colnames(coocCuesOutcomes)
    if (verbose) message("Completed calculations. Returning weight matrix.")
    flush.console()
    return(weightMatrix)
  }



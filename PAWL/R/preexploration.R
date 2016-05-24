preexplorationAMH <- function(target, nchains, niterations, proposal, verbose = TRUE){
    if (verbose) cat("Pre exploration to get log energy range\n")
    preexpparameters <- tuningparameters(nchains = nchains, niterations = niterations,
                                         adaptiveproposal = TRUE, storeall = FALSE)
    if (missing(proposal))
      preexp <- adaptiveMH(target, preexpparameters, verbose = verbose)
    else 
      preexp <- adaptiveMH(target, preexpparameters, proposal, verbose = verbose)
    burnin <- min(1000, niterations / 10)
    logtarget <- as.vector(preexp$alllogtarget[burnin:(niterations + 1),])
    LogEnergyRange <- range(-logtarget)
    LogEnergyQtile <- as.numeric(quantile(-logtarget, probs = 0.1))
    LogEnergyQtiles <- as.numeric(quantile(-logtarget, probs = c(0.1, 0.9)))
    return(list(LogEnergyRange = LogEnergyRange, 
                LogEnergyQtile = LogEnergyQtile,
                #SuggestedRange = c(LogEnergyQtile, LogEnergyRange[2]),
                SuggestedRange = c(LogEnergyQtiles[1], 
                                   LogEnergyQtiles[1] + 2 * (LogEnergyQtiles[2] - LogEnergyQtiles[1])),
                finalchains = preexp$finalchains))
}

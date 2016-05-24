# 2011-11-11 CJS Changed call from modelSetSeed to modelSetRN
# 2011-02-17 CJS InitialSeed in argument list changed to initial.seed
# 2011-02-17 CJS Removed reference to OpenBugs and WINBugs directory

run.openbugs <-
    function(modelFile,
             dataFile,
             dataList,
             initFiles,
             initVals,
             parameters,
             nChains,
             nIter,
             nBurnin,
             nSims,
             overRelax,
             initialSeed,
             working.directory,
             debug){

  ## The following call sequence is mainly based on the "openbugs"
  ## function in the R2WinBugs package.

  cat("\n\n*** Start of call to OpenBugs \n")
  cat("Working directory: ",working.directory,"\n")

  nIterPostBurnin <- nIter - nBurnin   # Number of iterations to run after burnin

  nThin <- round(nIterPostBurnin/nSims) # Thinning to obtain desired number of samples

  parametersToSave <- unique(c(parameters, "deviance")) # Always get the deviance

  cat("OpenBugs files created in:", working.directory, "\n")

  BRugs::modelCheck(modelFile)    # Check the model

  ## Format data
  bugsData <- BRugs::bugsData(dataList, fileName = dataFile, digits = 5)
  cat("Data files saved in ", dataFile, "\n")

  ## Write data file
  BRugs::modelData(dataFile)
  cat("Data loaded into model\n")

  ## Compile
  BRugs::modelCompile(nChains)
  cat("Model compiled\n")

  ## Set seed
  BRugs::modelSetRN(initialSeed)       ## BRUGS 0.7.1 which calls OPENBUS 3.2.1 Need to look into this
  ## BRugs::modelSetSeed(initialSeed)  ## BRugs 0.5-3.1
  cat("Initial seed for OpenBugs set to:", initialSeed, "\n")

  ## Format initial values and write to file
  inits <- BRugs::bugsInits(inits = initVals, numChains = nChains, fileName=initFiles)
  cat("Initial values generated in ",inits, "\n")

  ## Load initial values
  BRugs::modelInits(inits)
  cat("Initial values loaded.\n")

  ## Generate initial values for remaining parameters
  BRugs::modelGenInits()     # generate the initial values for any uninitialized variables
  cat("Initial values loaded into model\n")

  ## Generate burn-in samples
  cat("Burnin sampling has been started for ", nBurnin, " iterations.... \n")
  flush.console()
  for(iter in seq(1,nBurnin,round(nBurnin/20)+1)){  # generate a report about every 5% of the way
    cat('... Starting burnin iteration', iter,' which is about ',round(iter/nBurnin*100),
        "% of the burnin phase at ",date(),"\n")
    flush.console()
    BRugs::modelUpdate(round(nBurnin/20)+1, overRelax = overRelax)
  }
  cat("Burnin sampling completed \n")

  ## Initialize DIC
  BRugs::dicSet()      # turn on DIC computations
  on.exit(BRugs::dicClear(), add = TRUE)
  cat("DIC collection set \n")

  ## Set monitors
  BRugs::samplesSet(parametersToSave)
  cat("Nodes to monitor set\n")

  ## Generate remaining samples
  cat("Starting sampling after burnin for ", nChains," chain each with  a further ",
      nIterPostBurnin, " iterations. \n A thining rate of ", nThin,
      "will give about ", round(nIterPostBurnin/nThin), " posterior values in each chain... \n")

  for(iter in seq(1,nIterPostBurnin,round(nIterPostBurnin/20)+1)){
    cat('... Starting post-burnin iteration', iter,' which is about ',
        round(iter/nIterPostBurnin*100),"% of the post-burnin phase at ",date(),"\n")
    flush.console()
    BRugs::modelUpdate(round(nIterPostBurnin/nThin/20)+1, thin=nThin, overRelax = overRelax) # we do the thining on the fly
    BRugs::samplesCoda("*", stem=paste(working.directory,"/",sep=""))  # write out coda files for intermediate results
  }
  cat("Finished sampling after burnin and thinning \n")

  ## Extract the sampled values and create the bugs array
  cat("Extracting the sampled values\n")
  params <- BRugs::samplesMonitors("*")
  samples <- sapply(params,function(param){
      cat("     ",param,"\n")
      BRugs::samplesSample(param)
  })

  cat("Done\n")
  nSavedPerChain <- nrow(samples)/nChains
  samples.array <- array(samples, c(nSavedPerChain, nChains, ncol(samples)))
  dimnames(samples.array)[[3]] <- dimnames(samples)[[2]]
  DICOutput <- BRugs::dicStats()

  ## Save the results
  results<- as.bugs.array(sims.array = samples.array,
                          model.file = modelFile,
                          program = "OpenBUGS",
                          DIC = TRUE,
                          DICOutput = DICOutput,
                          n.iter = nIter,
                          n.burnin = nBurnin,
                          n.thin = nThin)

  results$Seed.initial <- NULL #initialSeed
  cat("Final dimension of saved simulation output is ", dim(results$sims.array), "\n")

                                        # Write CODA files to working directory
  BRugs::samplesCoda("*", stem=paste(working.directory,"/",sep=""))
  cat("Coda file created \n")


  ## Return results
  cat("\n\n*** Finished OpenBugs ***\n\n")
  results
}

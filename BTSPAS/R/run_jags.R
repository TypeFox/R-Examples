## 2014-09-01 CJS added code to dump out mcmc.list to coda files to match functionality of OpenBugs
##                Set the seed here.
## 2013-12-31 CJS added code to dump out data and initial values to files as in run.openbugs
## 2013-12-30 CJS changed program argument in as.bugs.array to JAGS
## 2013-09-22 sjb Created file. Copied from run_openbugs.R

run.jags <-
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
             debug=FALSE){

        cat("\n\n*** Start of call to JAGS \n")
        cat("Working directory: ",working.directory,"\n")

        nIterPostBurnin <- nIter - nBurnin   # Number of iterations to run after burnin
        
        nThin <- round(nIterPostBurnin/nSims) # Thinning to obtain desired number of samples
	#cat("In run_jags.R\n")
        #browser()

	## Set seed. We need to set a separate seed for each chain.
	## We start by setting the seed in R, and then generating nchain values between 1 and 10000000.
        cat("Initial seed for JAGS set to:", initialSeed, "\n")
	set.seed(initialSeed)
	initVals <- llply(initVals, function(x){
            # add to this list
            x$.RNG.seed <- round(runif(1, min=1, max=1000000))
	    x$.RNG.name <- "base::Wichmann-Hill"
	    cat("Random number seed for chain ", x$.RNG.seed, "\n")
	    x
        })

	# Dump out the data list (useful for debugging)
        #file.remove(dataFile)
	with(dataList, dump(names(dataList), file = dataFile))
	 
	#browser()

	# Dump out the initial values (useful for debugging)
        for (i in 1:nChains) {
            #file.remove(initFiles[i])
            initial.values <- initVals[[i]]
            with(initial.values, dump(names(initial.values), file = initFiles[i]))
	}
        
        #parametersToSave <- unique(c(parameters, "deviance")) # Always get the deviance
        parametersToSave <- unique(c(parameters, "deviance","pD"))

        ## Load the DIC model
        load.module("dic")
        
        ## Initialize JAGS model and run burnin
        jags <- jags.model(file=modelFile,
                           data=dataList,
                           inits=initVals,
                           n.chains=nChains,
                           n.adapt=nBurnin,
                           quiet=FALSE)

        ## Run sampling iterations
        samples <- coda.samples(model=jags,
                                variable.names=parametersToSave,
                                n.iter=nIterPostBurnin,
                                thin = nThin)

	## save the MCMC list as a coda file (to match OpenBugs)
	## taken from http://stackoverflow.com/questions/12078152/how-can-i-convert-an-mcmc-list-to-a-bugs-object
        s2 <- as.array(samples)
        lapply(seq_len(dim(s2)[3]), function(i) {
              write.table(cbind(rep(seq_len(nrow(s2[,,i])), ncol(s2)), c(s2[,,i])), 
                 paste0(working.directory, '/CODAchain', i, '.txt'),
                 row.names=FALSE, col.names=FALSE)
              })
        cat(paste(colnames(s2), 1+(seq_len(ncol(s2))-1) * nrow(s2), nrow(s2)*seq_len(ncol(s2))), 
               sep='\n', 
               file=file.path(working.directory, 'codaIndex.txt'))
        #browser()

        ## Convert output to array for input to as.bugs.array
        nSavedPerChain <- nrow(samples[[1]])
        samples.array.jags <- array(NA,c(nSavedPerChain, nChains, ncol(samples[[1]])))
        dimnames(samples.array.jags)[[3]] <- dimnames(samples[[1]])[[2]]
        
        for(i in 1:nChains){
#           there is subtle problem with two dimensional entries. OpenBUGS stores them in row major order; JAGS in column major
#           order and this screws up the entries. Convert the theta[] entries to the correct order by simple sort
#           before calling the as.bugs.array function 
	    Theta.cols <- grep("Theta", colnames(samples[[i]]))
	    if(length(Theta.cols)>0){
               temp        <-sapply(strsplit(colnames(samples[[i]])[Theta.cols],"\\["), "[", 2)
	       row.index   <-as.numeric(sapply(strsplit(temp                   ,","),"[",1))
	       temp        <-sapply(strsplit(temp,","), "[", 2)  # get the column index now
     	       col.index   <-as.numeric(sapply(strsplit(temp                   ,"\\]"),"[",1))
               samples[[i]][,Theta.cols] <- samples[[i]][,Theta.cols[order(row.index,col.index)]]
            } 
            samples.array.jags[,i,] <- samples[[i]]	    
        }

        ## Save the results
        results <- as.bugs.array(sims.array = samples.array.jags,
                                 model.file = "model.R",
                                 program = "JAGS",
                                 DIC = TRUE,
                                 ##DICOutput = DICOutput,
                                 n.iter = nIter,
                                 n.burnin = nBurnin,
                                 n.thin = nThin)
        ## Return results
        cat("\n\n*** Finished JAGS ***\n\n")

        return(results)
}

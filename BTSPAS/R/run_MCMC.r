# 2013-12-30 CJS - function to switch between the three samplers as needed. 
#                  This way we wont't have to modify much code (hopefully)

run.MCMC <-
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
             overRelax=FALSE,
             initialSeed,
             working.directory,
             debug=FALSE,
	     engine=c("jags","openbugs")[1]){

if(tolower(engine) %in% c("openbugs")){
    results <- run.openbugs(modelFile=modelFile,
                            dataFile=dataFile,
                            dataList=dataList,
                            initFiles=initFiles,
                            initVals=initVals,
                            parameters=parameters,
                            nChains=nChains,
                            nIter=nIter,
                            nBurnin=nBurnin,
                            nSims=nSims,
                            overRelax=overRelax,
                            initialSeed=initialSeed,
                            working.directory=working.directory,
                            debug=debug)
   } 
else if(tolower(engine) == "jags"){
        results <- run.jags(modelFile=modelFile,
                           dataFile=dataFile,
                           dataList=dataList,
                           initFiles=initFiles,
                           initVals=initVals,
                           parameters=parameters,
                           nChains=nChains,
                           nIter=nIter,
                           nBurnin=nBurnin,
                           nSims=nSims,
                           overRelax=overRelax,
                           initialSeed=initialSeed,
                           working.directory=working.directory,
                           debug=debug)
    } 
else{ stop(paste("Unknown BUGS engine:",engine,".\n")) }

return(results)
} # end of function


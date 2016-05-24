###########################################################################################
# This file contains the functions needed to check the MCMC chains for convergence,
#   length of run time, etc.
###########################################################################################

checkRunTime = function(maxIter, systemtime){
    
    currenttime = unclass(Sys.time())
    estruntime = round((currenttime-unclass(systemtime)-.05)*min(c(maxIter, 100000))/100)
    units = 'seconds'
    if(estruntime > 60){
        units = 'minutes'
        estruntime = round(estruntime/60)
	}
    if(estruntime > 60){
        units = 'hours'
        estruntime = ceiling(estruntime/60)
	}
	writeLines(paste('Estimated total run time is', estruntime, units, '\n'))
	writeLines('To shorten the run time, re-run with fewer iterations or a smaller number of bins. \n')
}

checkMCMCParams = function(best.thin.found, assessIndex, thin, burnIn, round.values, outfilename,
convergence, fix.burnIn){

    testData = read.table(paste(outfilename, '/MCMCchains.txt', sep = ''), header = TRUE)

    # Test autocorrelation and find best thinning value
	while(best.thin.found == FALSE){
        newthin = thinAmt(thin, mcmc(testData[,assessIndex], thin = thin, start = burnIn))
        if(newthin == thin){    best.thin.found = TRUE
        } else {
            testData = testData[seq(1, nrow(testData), by = newthin/thin),]
            thin = newthin
		}
    }
	keepnames = colnames(testData)
	testData = cbind(testData[,1], sapply(2:ncol(testData), function(x) round(testData[,x], round.values[x-1])))
	colnames(testData) = keepnames
	write.table(testData, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
			
	# Test convergence
	enough.data.left = TRUE
	if(nrow(testData) < 500){ enough.data.left = FALSE	}
	newburnIn = burnIn
    while(convergence == FALSE & enough.data.left == TRUE){
        convergence = convergeFxn(mcmc(testData[,assessIndex], thin = thin, start = burnIn))
		if(fix.burnIn == FALSE & convergence == FALSE & nrow(testData) > 30000/thin){
            testData = testData[-c(1:(20000/thin)),]
            newburnIn = newburnIn + 20000
		} else { enough.data.left = FALSE	}
    }
	if(convergence == TRUE & fix.burnIn == TRUE){
        burnIn = newburnIn
		keepnames = colnames(testData)
		testData = cbind(testData[,1], sapply(2:ncol(testData), function(x) round(testData[,x], round.values[x-1])))
		colnames(testData) = keepnames
		write.table(testData, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
    }
    returndata = list(thin = thin, burnIn = burnIn, convergence = convergence,
                        best.thin.found = best.thin.found)
}


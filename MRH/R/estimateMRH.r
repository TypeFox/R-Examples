###########################################################################################
# This is the MRH wrapper
###########################################################################################

estimateMRH = function(formula, data, M, numberBins, maxStudyTime, outfolder = 'MRHresults', 
prune = FALSE, prune.alpha = 0.05, prune.levels = NULL, burnIn = 50000, maxIter = 500000, 
thin = 10, Rmp.init = NULL, a.init = 10, lambda.init, beta.init = NULL, k.fixed,
gamma.fixed, GR = FALSE, convGraphs = TRUE, fix.burnIn = FALSE, fix.thin = FALSE, fix.max = FALSE,
continue.chain = FALSE){

	options(digits = 16)
	systemtime = Sys.time()
	sysseconds = format(systemtime, "%H%M%S")
	dir.create(outfolder, showWarnings = FALSE)  
	
	writeLines(paste('\n', 'MCMC routine running.  Calculating estimated runtime for', 
					 min(100000, maxIter), 'iterations...', '\n'))
	
	nonph.model = FALSE
	if(!missing(data)){		model = model.frame(formula, data = data)
	} else {	model = model.frame(formula)	}
	if(!inherits(model.extract(model, 'response'), "Surv")){	stop("Response must be a survival object")	}
	specialcols = attr(terms(formula, specials = "nph"), "specials")$nph
	if(length(specialcols) > 0){ nonph.model = TRUE	}
	Responses = as.matrix(model.extract(model, 'response'))
	Ti = Responses[,1]
	delta = Responses[,2]
	if(nonph.model == TRUE){
		tempmodel = model.frame(model)
		if(ncol(tempmodel) > 2){
			Xmodel = model.frame(model[,-specialcols])
			Xmatrix = model.matrix(attr(Xmodel, 'terms'), data = Xmodel)
			namesXmat = colnames(Xmatrix)
			Xmatrix = as.matrix(Xmatrix[,-1])
			colnames(Xmatrix) = namesXmat[-1]
		} else {	Xmatrix = NULL	}
		numNPHcovs = 1
		if(length(specialcols) == 1){
			X.nonproportional = as.matrix(model[,specialcols])
			colnames(X.nonproportional) = strsplit(unlist(strsplit(names(model)[specialcols], "\\("))[2], "\\)")
			numPHcovs = ncol(Xmatrix)
		} else {
			tempmat = model.matrix(~0 + as.factor(model[,specialcols[1]]):as.factor(model[,specialcols[2]]))
			newfactorvar = as.data.frame(cbind(1:nrow(tempmat), 
											   matrix(rowSums(tempmat*matrix(1:ncol(tempmat), byrow = TRUE, nrow = nrow(tempmat), 
																			 ncol = ncol(tempmat))), ncol = 1)))
			names(newfactorvar) = c('orderval', 'factorlevel')
			firstvals = unique(model[,specialcols[1]])
			firstvals = firstvals[order(firstvals)]
			secvals = unique(model[,specialcols[2]])
			secvals = secvals[order(secvals)]
			colctr = 1
			names.format = as.data.frame(cbind(1:ncol(tempmat), NA))
			names(names.format) = c('factorlevel', 'factorname')
			for(colctr2 in 1:length(secvals)){ for(colctr1 in 1:length(firstvals)){	
				names.format$factorname[colctr] = paste(strsplit(unlist(strsplit(names(model)[specialcols[1]], "\\("))[2], "\\)"), firstvals[colctr1], '.', 
														strsplit(unlist(strsplit(names(model)[specialcols[2]], "\\("))[2], "\\)"), secvals[colctr2], sep = '')
				colctr = colctr + 1
			}}
			nph.factors = merge(newfactorvar, names.format, by = 'factorlevel')
			nph.factors = nph.factors[order(nph.factors$orderval),]
			if(length(specialcols) > 2){  
				for(i in 3:length(specialcols)){
					tempmat = model.matrix(~0 + as.factor(nph.factors$factorlevel):as.factor(model[,specialcols[i]]))
					newfactorvar = as.data.frame(cbind(1:nrow(tempmat), 
													   matrix(rowSums(tempmat*matrix(1:ncol(tempmat), byrow = TRUE, nrow = nrow(tempmat), 
																					 ncol = ncol(tempmat))), ncol = 1)))
					names(newfactorvar) = c('orderval', 'factorlevel')
					firstvals = unique(nph.factors$factorname)
					firstvals = firstvals[order(unique(nph.factors$factorlevel))]
					secvals = unique(model[,specialcols[i]])
					secvals = secvals[order(secvals)]
					colctr = 1
					names.format = as.data.frame(cbind(1:ncol(tempmat), NA))
					names(names.format) = c('factorlevel', 'factorname')
					for(colctr2 in 1:length(secvals)){ for(colctr1 in 1:length(firstvals)){	
						names.format$factorname[colctr] = paste(firstvals[colctr1], '.', 
																strsplit(unlist(strsplit(names(model)[specialcols[i]], "\\("))[2], "\\)"), secvals[colctr2], sep = '')
						colctr = colctr + 1
					}}
					nph.factors = merge(newfactorvar, names.format, by = 'factorlevel')
					nph.factors = nph.factors[order(nph.factors$orderval),]
				}
			}
			X.nonproportional = matrix(nph.factors[,-(1:2)], ncol = 1)
			colnames(X.nonproportional) = paste(strsplit(unlist(strsplit(names(model)[specialcols], "\\("))[2], "\\)"), sep = '.')
		}
	} else {
		Xmatrix = model.matrix(attr(model, 'terms'), data = model)
		if(ncol(Xmatrix) == 1){	Xmatrix = NULL	
		} else { 
			namesXmat = colnames(Xmatrix)
			Xmatrix = as.matrix(Xmatrix[,-1])	
			if(ncol(Xmatrix) == 1){	colnames(Xmatrix) = namesXmat[2]	}
		}
		numPHcovs = ncol(Xmatrix)
		numNPHcovs = 0
		numHazards = 1
	}
	## Perfect the names ##
	if(!is.null(Xmatrix)){	colnames(Xmatrix) = gsub("`", "", colnames(Xmatrix))	}
	if(nonph.model == TRUE){	colnames(X.nonproportional) = gsub("`", "", colnames(X.nonproportional))	}

	###########################################################################################
	#		Check user-entered information
	###########################################################################################	
	# Check that either the number of bins or M has been entered, and if the number
	# of bins is entered, it is a power of two.
	if(missing(M) & missing(numberBins)){	stop("Must enter number of bins or M")	}
	if(!missing(numberBins)){ 
		M = log(numberBins)/log(2)
		if(M %% 1 != 0){ stop("Number of bins must be a factor of 2")	}
	}
	
	# Check bounds on user-entered initialized parameters
	if(length(which(0 > Rmp.init | Rmp.init > 1)) > 0){ stop("All Rmp values must be between 0 and 1.") }
	if(!is.null(Rmp.init) & length(Rmp.init) != 2^M-1){ stop("Length of initialized Rmp vector should equal 2^M-1.") }
	if(a.init < 0){ stop("Initial value for 'a' must be greater than 0") }
	if(!missing(lambda.init)){	if(lambda.init < 0){ stop("Initial value for lambda must be greater than 0") }}
	# Check that the pruning vector or matrix has the correct length/dimensions if it is entered.
	if(!is.logical(prune)){ 
		if(nonph.model == FALSE & length(prune) != (2^M-1)){	
			stop("Pruning indicator (prune) is incorrect.  Each Rmp must have a corresponding 1 or 0.")	
		} else if (nonph.model == TRUE){
			if((ncol(prune) > 1 & (ncol(prune) != length(names(table(X.nonproportional))) | nrow(prune) != 2^M-1)) | 
			   (is.null(ncol(prune)) & length(prune) != 2^M-1)){
					stop("Pruning indicator is incorrect.  Indicator should either be a vector used for
						 all hazard rates for all groups, or a matrix with a pruning indicator in each column for each group.")
			}
		}
	}
	if(!is.null(beta.init)){ if(length(beta.init) != numPHcovs){ 
		stop("Number of initialized beta values does not match the number of covariates in the model.") 
	}}

	##########################################################################
	# If user has selected the Gelman-Rubin option, then set the burn-in, 
	# max iterations, and thinning values.  If the maximum number of iterations
	# is less than 100,000, then change the convergence check for that. 
	# If the maximum number of iterations is less than 10,000, make sure to 
	# write out data set before that.
	##########################################################################
	if(GR == TRUE){		fix.burnIn = fix.max = fix.thin = TRUE	}
	writenum = 10000
	checknum = 100000
	if(maxIter < 100000){	
		checknum = maxIter
		if(maxIter < 10000){	writenum = maxIter	} 
	}
	if(burnIn >= maxIter){	burnIn = round(.5*maxIter)	}
	
	###########################################################################################
	#		Calculate the censoring time, set pruning indicator, etc.
	###########################################################################################
	if(missing(maxStudyTime)){	censortime = max(Ti)
	} else {
		censortime = maxStudyTime
		delta[delta > maxStudyTime] = 0
		Ti[Ti > maxStudyTime] = maxStudyTime
	}
	if(nonph.model == TRUE){
		IDtypes = names(table(X.nonproportional))
		IDtypes = IDtypes[order(IDtypes)]
		sepHazNames = IDtypes
		numHazards = length(IDtypes)
		indices = list(which(X.nonproportional == IDtypes[1]))
		for(hazCtr in 2:numHazards){	indices = c(indices, list(which(X.nonproportional == IDtypes[hazCtr])))	}
	}
		
	##########################################################################
	# If user has selected the Gelman-Rubin option, then set the burn-in, 
	# max iterations, and thinning values.  
	##########################################################################
	if(GR == TRUE){		fix.burnIn = fix.max = fix.thin = TRUE	}

    ################################################################################
    # If user wants to continue the previous chain, then set the initial values
    # from the output folder, and use the previous burn-in and thinning values as
    # the initial starting parameter values.
    #################################################################################
    loopctr.start = 1
    if(continue.chain == TRUE){
        mcmcinfo = read.table(paste(outfolder, '/MCMCchains.txt', sep = ''), header = TRUE)
        # Need to remove the last column that contains the log-likelihood calculations
        write.table(mcmcinfo[,-ncol(mcmcinfo)], paste(outfolder, '/MCMCchains.txt', sep = ''),
                    row.names = FALSE)
        burnIn = 0
        thin = mcmcinfo[2,1]-mcmcinfo[1,1]
        maxIter = mcmcinfo[nrow(mcmcinfo),1]+maxIter
        loopctr.start = mcmcinfo[nrow(mcmcinfo),1]+1
        
        lastline = matrix(as.numeric(mcmcinfo[nrow(mcmcinfo),]), nrow = 1)
        colnames(lastline) = names(mcmcinfo)
        
        betavals = which(substr(colnames(lastline), 1, 4) == 'beta')
        if(numHazards > 1){
            badbetavals = which(grepl('.bin1', colnames(lastline)) == 'TRUE')[1]
            if(length(badbetavals) > 0 & length(betavals) > 0){
                if(betavals[1] == badbetavals){ betavals = NULL
                } else {    betavals = betavals[1:which(betavals == (badbetavals-1))] }
            }
        }
        if(length(betavals) > 0){   beta.init = lastline[,betavals]  }
        
        rmpvals = which(substr(colnames(lastline), 1, 3) == 'Rmp')
        Rmp.init = matrix(lastline[,rmpvals], ncol = numHazards)
        avals = which(substr(colnames(lastline), 1, 1) == 'a')
        a.init = lastline[,avals]
        
        lambdavals = which(substr(colnames(lastline), 1, 6) == 'lambda')
        lambda.init = lastline[,lambdavals]
        
        gammavals = which(substr(colnames(lastline), 1, 5) == 'gamma')
        if(length(gammavals) > 0){
            gamma.init = matrix(lastline[,gammavals], ncol = numHazards)
            gamma.fixed = FALSE
        }
        
        kvals = which(substr(colnames(lastline), 1, 1) == 'k')
        if(length(kvals) > 0){
            k.init = lastline[,kvals]
            k.fixed = FALSE
        }
    }

	###########################################################################################
	#		Create the pruning variable based on user entered information
	###########################################################################################
	prune.indicator = NULL
	if(is.logical(prune) & prune == TRUE){
		if(M < 2){ stop("Pruning not available for less than 2 bins.")		}
		# If it is not the non-proportional model, then check for errors in the prune levels, 
		# and set them appropriately
		if(nonph.model == FALSE){
			if(!is.null(prune.levels) & length(prune.levels) > 1){	
				stop("prune.level value is incorrect. Only enter 1 value for pruning level")
			} else if (is.null(prune.levels)) {	prune.levels = M	}
			prune.indicator = Prune(Ti, delta, M, prune.levels, censortime, alpha = prune.alpha)
		} else {
			if(!is.null(prune.levels) & length(prune.levels) > 1 & length(prune.levels) != numHazards){
					stop("prune.level value is incorrect. 
						 Enter 1 value for pruning all hazards, or enter a value for each hazard")
			}
			if(is.null(prune.levels)){ prune.levels = rep(M, numHazards)	
			} else if(!is.null(prune.levels) & length(prune.levels) == 1){
				prune.levels = rep(prune.levels, numHazards)
			}
			prune.indicator = Prune(Ti[indices[[1]]], delta[indices[[1]]], M, prune.levels[1],
										censortime, alpha = prune.alpha)
			for(hazCtr in 2:numHazards){ 
				prune.indicator = cbind(prune.indicator, 
										Prune(Ti[indices[[hazCtr]]], delta[indices[[hazCtr]]], M, 
											  prune.levels[hazCtr], censortime, alpha = prune.alpha))
			}
		}
	} else if(!is.logical(prune)) {	prune.indicator = prune	}

	###########################################################################################
	#		Set the k and gamma values, checking for user errors.
	###########################################################################################
	if(missing(lambda.init)){	lambda.init = NULL	}
	
	# If the user does not enter anything for the k.fixedor if it is set to TRUE, set k.fixed = 0.5
	if(missing(k.fixed)){	k.fixed = rep(0.5, numHazards)	
	} else if (k.fixed == TRUE){	k.fixed = rep(0.5, numHazards)
	# If the user only enters one k.fixed value for the NPH case, set that value to all hazards
	} else if (!is.logical(k.fixed) & length(k.fixed) == 1 & numHazards > 1){ k.fixed = rep(k.fixed, numHazards)
	# If more than one k is entered, but it does not match the number of hazards, send error message
	} else if (!is.logical(k.fixed) & length(k.fixed) > 1 & length(k.fixed) != numHazards){ 
			stop("Enter one value of k for all hazards or a vector of k values with one for each hazard")
	# If the user enters a value of k less than zero, send error message.
	} else if(!is.logical(k.fixed) & length(which(k.fixed <= 0)) > 0){	stop("Values of k must be greater than 0.")	}
	
	# If the user does not enter a value for gamma or if it is set to TRUE, set gamma.fixed at 0.5
	if(missing(gamma.fixed)){	gamma.fixed = matrix(0.5, nrow = 2^M-1, ncol = numHazards)
	} else if (gamma.fixed[1] == TRUE){	gamma.fixed = matrix(0.5, nrow = 2^M-1, ncol = numHazards)
	# If the user only enters one vector of gammas (for the NPH case), set that vector to all hazards		
	} else if (!is.logical(gamma.fixed[1]) & length(gamma.fixed) == 2^M-1 & numHazards > 1){ 
		gamma.fixed = matrix(gamma.fixed, nrow = 2^M-1, ncol = numHazards)
	# If the number of gammas entered does not match the number of hazards, sent and error messge
	} else if (!is.logical(gamma.fixed[1]) & length(gamma.fixed) != (2^M-1)*numHazards){
		stop(paste("Incorrect number of gamma values entered:", 2^M-1, "fixed gamma values needed"))
	} else if (!is.logical(gamma.fixed[1]) & length(which(gamma.fixed < 0 | gamma.fixed > 1)) > 0){
		stop("All values of gamma must be between 0 and 1")	}
	
	###########################################################################################
	#		Call the correct routine
	###########################################################################################
	#################################### PH MODEL ####################################
	if(nonph.model == FALSE){
		if(k.fixed == FALSE & gamma.fixed[1] == FALSE){ sampletype = 'kandg'
        } else if(k.fixed == FALSE & gamma.fixed[1] != FALSE){  sampletype = 'konly'
        } else if(k.fixed != FALSE & gamma.fixed[1] == FALSE){  sampletype = 'gonly'
        } else {    sampletype = 'simple'   }
        outputdata = PHshell(Mval = M, Ti = Ti, delta = delta, X = Xmatrix, outfilename = outfolder,
                                censortime = censortime, prune.indc = prune.indicator, burnIn = burnIn,
                                maxIter = maxIter, thin = thin, RmpInit = Rmp.init, a.init = a.init,
                                lambda.init = lambda.init, beta.init = beta.init, k = k.fixed,
                                gamma.mp = gamma.fixed, GR = GR, fix.burnIn = fix.burnIn,
                                fix.thin = fix.thin, fix.max = fix.max, writenum = writenum,
                                sysseconds = sysseconds, systemtime = systemtime, checknum = checknum,
                                sampletype = sampletype, continue.chain = continue.chain, k.init = k.init,
                                gamma.init = gamma.init, loopctr.start = loopctr.start)
                                
		######## Calculate the log-likelihood values, the information criteria values, make the convergence graphs, make the summaries
		finaldata = read.table(paste(outfolder, '/MCMCchains.txt', sep = ''), header = TRUE)
		# Log likelihood and AIC, BIC, DIC
		ICinfo = loglikePH(finaldata = finaldata, numEstParams = outputdata$numberofEstParams, censortime = censortime, Mval = M,
                            numParams = outputdata$numParams, delta = delta, failBin = outputdata$failBin,
							inBin = outputdata$inBin, outfilename = outfolder, X = Xmatrix)
        nameshazgroups = NULL
        numHazards = 1
	#################################### NPH MODEL ####################################
	# Non-proportional MRH
	} else {
		if(k.fixed[1] == FALSE & gamma.fixed[1] == FALSE){ sampletype = 'kandg'
        } else if(k.fixed[1] == FALSE & gamma.fixed[1] != FALSE){ sampletype = 'konly'
        } else if(k.fixed[1] != FALSE & gamma.fixed[1] == FALSE){ sampletype = 'gonly'
        } else {    sampletype = 'simple'   }
		outputdata = NPHBAshell(Mval = M, Ti = Ti, Xfixed = Xmatrix, XNPH = X.nonproportional,
								 prune.indc = prune.indicator, delta = delta, outfilename = outfolder, 
								 burnIn = burnIn, maxIter = maxIter, thin = thin, RmpInit = Rmp.init,
								 a.init = a.init, lambda.init = lambda.init, beta.init = beta.init, k = k.fixed, gamma.mp = gamma.fixed, 
								 censortime = censortime, GR = GR, convGraphs = convGraphs, fix.burnIn = fix.burnIn,
								 fix.thin = fix.thin, fix.max = fix.max, writenum = writenum, 
								 sysseconds = sysseconds, systemtime = systemtime, checknum = checknum,
                                 sampletype = sampletype, continue.chain = continue.chain, k.init = k.init,
                                 gamma.init = gamma.init, loopctr.start = loopctr.start)

		######## Calculate the log-likelihood values, the information criteria values, make the convergence graphs, make the summaries
        IDtypes = names(table(X.nonproportional))
        IDtypes = IDtypes[order(IDtypes)]
        numHazards = length(IDtypes)
        indices = list(which(X.nonproportional == IDtypes[1]))
        for(hazCtr in 2:numHazards){
            indices = c(indices, list(which(X.nonproportional == IDtypes[hazCtr])))
        }
		finaldata = read.table(paste(outfolder, '/MCMCchains.txt', sep = ''), header = TRUE)
		ICinfo = loglikeNPH(finaldata, outputdata$numberofEstParams, censortime, M, outputdata$numParams, delta,
                            outputdata$failBin, outputdata$inBin, outfolder, numHazards = numHazards,
                            Xfixed = Xmatrix, numPHParams = outputdata$numPHParams, indices = indices)
        nameshazgroups = ICinfo$namesHazGroups
	}

	########################### Plot the MCMC chain diagnostics ####################
	if(convGraphs == TRUE){
		convergenceGraphs(finaldata[,outputdata$assessIndex], paste(outfolder, '/convergenceGraphs', sep = ''), burnIn, thin)
	}
	########################### Summaries of MCMCchains	###########################
	summarydata = AnalyzeComplete(M = M, censortime = censortime, finaldataset = finaldata,
                                    outfilename = paste(outfolder, '/', sep = ''), numhazards = numHazards,
                                    namesHazGroups = nameshazgroups)
	if(nonph.model == TRUE){
		if(outputdata$numPHParams > 0){
			GraphNPbetas(M, dests = summarydata$d, np.betaests = summarydata$beta[-(1:outputdata$numPHParams),], 
						 numhazgroups = numHazards, hazgroupnames = outputdata$namesHazGroups, censortime = censortime, 
						 file = paste(outfolder, '/', sep = ''))
		} else {
			GraphNPbetas(M , dests = summarydata$d, np.betaests = summarydata$beta, 
						 numhazgroups = numHazards, hazgroupnames = outputdata$namesHazGroups, censortime = censortime, 
						 file = paste(outfolder, '/', sep = ''))
		}
	}		
	
	if(outputdata$convergence == FALSE){
		warning(paste("Algorithm has not yet converged after ", outputdata$TotalIters, " MCMC iterations. 
					  Parameter estimates may not be reliable. ", sep = ''))
	}
	if(GR == FALSE){ gr.used = 'Option not used' } else { gr.used = 'Option Used'	}
    # Calculate final burn-in and thin values used
    burnIn.final = finaldata[1,1]
    thin.final = finaldata[2,1]-finaldata[1,1]
	returndata = c(summarydata, ICinfo, list(burnIn = burnIn.final, thin = thin.final, TotalIters = outputdata$loopctr-1,
				   convergence = outputdata$convergence, gelman.rubin.used = gr.used, fix.thin = fix.thin,
				   fix.burnIn = fix.burnIn, fix.max = fix.max, initialValues = outputdata$initialValues,
				   k.fixed = k.fixed, gamma.fixed = gamma.fixed, 
				   runtime = paste(round((unclass(Sys.time()) - unclass(systemtime))/(60*60), 2), 'hours'), 
				   outfolder = outfolder, maxStudyTime = censortime))
	class(returndata) = "MRH"
	return(returndata)
		
}

######################################################################################
# Contains extra functions beyond MRH estimation that the user may want for analyzing
#	results.
######################################################################################

##################################################################################
# binOptimize: This function assists the user in finding the optimal number of bins
#	and the optimal time unit for analysis with the MRH library, and returns the
#	new times and censoring variable.
##################################################################################
FindBinWidth = function(time, delta, time.unit, maxStudyTime){

	if(missing(maxStudyTime)){	
		censortime = max(time)	
	} else {
		censortime = maxStudyTime
		print(noquote(paste('With the given study length, ', length(which(time > censortime & delta == 1)),
					' extra observed failures will be censored.')))
		cat(' \n')
		time[time > censortime] = censortime
	}

	if(time.unit == 's'){	
		wordtime.unit = 'seconds'
		censortime = censortime
	} else if(time.unit == 'min'){
		wordtime.unit = 'minutes'
		censortime = censortime*60
	} else if(time.unit == 'h'){
		wordtime.unit = 'hours'
		censortime = censortime*60*60
	} else if(time.unit == 'd'){
		wordtime.unit = 'days'
		censortime = censortime*24*60*60
	} else if(time.unit == 'w'){	
		wordtime.unit = 'weeks'
		censortime = censortime*7*24*60*60
	} else if(time.unit == 'mon'){	
		wordtime.unit = 'months'
		censortime = censortime*365.25/12*24*60*60
	} else if(time.unit == 'y'){	
		wordtime.unit = 'years'
		censortime = censortime*365.25*24*60*60
	}	
		
	print(noquote(paste('The mean failure time is ', round(mean(time), 3), wordtime.unit)))
	print(noquote(paste('The median failure time is ', round(median(time), 3), wordtime.unit)))
	print(noquote(paste('The range of failure times is ', round(min(time), 3), 'to',
				round(max(time), 3), wordtime.unit)))
	cat(' \n')
	
	Mpossible = 2:10
	time.perbin = censortime/2^Mpossible
	
	orig.perbin = as.data.frame(cbind(time.perbin, time.perbin/60, time.perbin/(60*60),
									  time.perbin/(60*60*24), time.perbin/(60*60*24*7),
									  time.perbin/(60*60*24*365.25/12), time.perbin/(60*60*24*365.25)))
	names(orig.perbin) = c('secs', 'mins', 'hours', 'days', 'weeks', 'months', 'years')
	row.names(orig.perbin) = paste('M', Mpossible, sep = '')
	
	return(orig.perbin)
	
}

###########################################################################################
# The AnalyzeMultiple function calculates the estimates and bounds 
###########################################################################################
AnalyzeMultiple = function(datalist, fileNames, alpha.level = 0.05, maxStudyTime, GR = TRUE){
	
	censortime = maxStudyTime
	
	if(missing(datalist) & missing(fileNames)){ stop("Need a list of data sets or the file names 
													 for the data sets.") }
	if(missing(fileNames)){	numsets = length(datalist)	
	} else { numsets = length(fileNames)	}

	if(missing(fileNames)){	pdata = MRH(datalist[[1]])
	} else {	pdata = MRH(read.table(fileNames[1], header = TRUE))	}

	nonph.model = TRUE
	if(substr(colnames(pdata)[2], 3, 3) == ""){	nonph.model = FALSE	}
	pdata = pdata[,-c(1,ncol(pdata))]
	if(nonph.model == FALSE){		
		M = log(length(which(substr(colnames(pdata), 1, 1) == 'd')))/log(2)
		numhazards = 1
		numphbetas = ncol(pdata)-2^M*2
		if(numphbetas > 0){ phbetanames = names(pdata)[(2^M+1):(ncol(pdata)-2^M-1)]	}
		length.uniques = sapply(1:ncol(pdata), function(x) length(unique(pdata[,x])))
		pruned.Rmps = which(length.uniques == 1 & pdata[1,] == 0.5 & substr(colnames(pdata), 1, 3) == 'Rmp')
		if(length(pruned.Rmps) > 0){
            if(length(which(substr(colnames(pdata), 1, 3) == 'gam')) > 0){
                pruned.Rmps = c(pruned.Rmps, pruned.Rmps+2^M-1)
            }
			gr.keep.index = c((2^M+1):ncol(pdata))[-c(pruned.Rmps-2^M)]
		} else {	gr.keep.index = c((2^M+1):ncol(pdata))	}
	} else {
		firstgroup = substr(colnames(pdata)[2], 4, 4)
		dnames = colnames(pdata)[which(substr(colnames(pdata), 1, 1) == 'd')]
		numbins = length(which(unlist(strsplit(dnames, '_'))[1:length(dnames)*2] == firstgroup))
		M = log(numbins)/log(2)
		hazgroups = unique(unlist(strsplit(dnames, '_'))[1:length(dnames)*2])
		numhazards = length(hazgroups)
		numphbetas = ncol(pdata)-(2^M*2*numhazards + 2^M*(numhazards-1))-1
		if(numphbetas > 0){ phbetanames = names(pdata)[(2^M*numhazards+1):(ncol(pdata)-2^M*numhazards-1)]	}
		length.uniques = sapply(1:ncol(pdata), function(x) length(unique(pdata[,x])))
		pruned.Rmps = which(length.uniques == 1 & pdata[1,] == 0.5 & substr(colnames(pdata), 1, 3) == 'Rmp')
		if(length(pruned.Rmps) > 0){
            if(length(which(substr(colnames(pdata), 1, 3) == 'gam')) > 0){
                pruned.Rmps = c(pruned.Rmps, pruned.Rmps+(2^M-1)*numhazards)
            }
			gr.keep.index = c((2^M*numhazards+1):ncol(pdata))[-c(pruned.Rmps-2^M*numhazards)]
			keep.names = colnames(pdata)[gr.keep.index]
			beta.nph.start = which(substr(keep.names, nchar(keep.names)-3, nchar(keep.names))=='bin1')[1]
			gr.keep.index = gr.keep.index[-c(beta.nph.start:(beta.nph.start-1+2^M*(numhazards-1)))]
		} else {
			gr.keep.index = c((2^M*numhazards+1):ncol(pdata))
			keep.names = colnames(pdata)[gr.keep.index]
			beta.nph.start = which(substr(keep.names, nchar(keep.names)-3, nchar(keep.names))=='bin1')[1]
			gr.keep.index = gr.keep.index[-(beta.nph.start:(beta.nph.start-1+2^M*(numhazards-1)))]
		}
	}
	
	hazrate = survrate = cumulhaz = d = H = Rmp = beta = betanph = k = gamma = mcmctemplist = burnIn = numIters = NULL
	gr.datalist = NULL
	for(filectr in 1:numsets){
		
		#### Store the parameter file in pdata ###
		if(missing(fileNames)){		temp.pdata = datalist[[filectr]]
		} else {	temp.pdata = read.table(fileNames[filectr], header = TRUE)	}
		pdata = MRH(temp.pdata)
		gr.datalist = c(gr.datalist, list(temp.pdata[,-c(1,ncol(temp.pdata))][,gr.keep.index]))
	
		### Keep track of the burn in for each data set.	
		### This is needed for the Gelman-Rubin test.
		burnIn = c(burnIn, pdata[1,1])
		numIters = c(numIters, pdata[nrow(pdata),1])
		thinval = pdata[2,1]-pdata[1,1]
		
		#### Calculate the medians and alpha-level quantiles ####	
		pdata.summary = summary(pdata, alpha.level = alpha.level, maxStudyTime = maxStudyTime)
		hazrate = cbind(hazrate, as.matrix(pdata.summary$hazardRate))
		survrate = cbind(survrate, as.matrix(pdata.summary$SurvivalCurve))
		cumulhaz = cbind(cumulhaz, as.matrix(pdata.summary$CumulativeHazard))
		if('beta' %in% names(pdata.summary)){	beta = cbind(beta, as.matrix(pdata.summary$beta))	
		} else if('betaPH' %in% names(pdata.summary)){
			beta = cbind(beta, as.matrix(pdata.summary$betaPH))
			betanph = cbind(betanph, as.matrix(pdata.summary$betaNPH))
		}
		d = cbind(d, as.matrix(pdata.summary$d))
		H = cbind(H, as.matrix(pdata.summary$H))
		Rmp = cbind(Rmp, as.matrix(pdata.summary$Rmp))
		if('k' %in% names(pdata.summary)){	k = cbind(k, as.matrix(pdata.summary$k))	}
		if('gamma' %in% names(pdata.summary)){	gamma = cbind(gamma, as.matrix(pdata.summary$gamma))	}
			
	}
	
	hazrate.full = cbind(sapply(1:nrow(hazrate), function(x) median(hazrate[x,(1:numsets-1)*3+1])),
						 sapply(1:nrow(hazrate), function(x) median(hazrate[x,(1:numsets-1)*3+2])),
						 sapply(1:nrow(hazrate), function(x) median(hazrate[x,(1:numsets-1)*3+3])))
	survrate.full = cbind(sapply(1:nrow(survrate), function(x) median(survrate[x,(1:numsets-1)*3+1])),
						  sapply(1:nrow(survrate), function(x) median(survrate[x,(1:numsets-1)*3+2])),
						  sapply(1:nrow(survrate), function(x) median(survrate[x,(1:numsets-1)*3+3])))
	cumulhaz.full = cbind(sapply(1:nrow(cumulhaz), function(x) median(cumulhaz[x,(1:numsets-1)*3+1])),
						  sapply(1:nrow(cumulhaz), function(x) median(cumulhaz[x,(1:numsets-1)*3+2])),
						  sapply(1:nrow(cumulhaz), function(x) median(cumulhaz[x,(1:numsets-1)*3+3])))
	
	d.full = cbind(sapply(1:nrow(d), function(x) median(d[x,(1:numsets-1)*3+1])),
				   sapply(1:nrow(d), function(x) median(d[x,(1:numsets-1)*3+2])),
				   sapply(1:nrow(d), function(x) median(d[x,(1:numsets-1)*3+3])))
	H.full = cbind(sapply(1:nrow(H), function(x) median(H[x,(1:numsets-1)*3+1])),
				   sapply(1:nrow(H), function(x) median(H[x,(1:numsets-1)*3+2])),
				   sapply(1:nrow(H), function(x) median(H[x,(1:numsets-1)*3+3])))
	Rmp.full = cbind(sapply(1:nrow(Rmp), function(x) median(Rmp[x,(1:numsets-1)*3+1])),
					 sapply(1:nrow(Rmp), function(x) median(Rmp[x,(1:numsets-1)*3+2])),
					 sapply(1:nrow(Rmp), function(x) median(Rmp[x,(1:numsets-1)*3+3])))
	if('beta' %in% names(pdata.summary)){
		beta.full = cbind(sapply(1:nrow(beta), function(x) median(beta[x,(1:numsets-1)*3+1])),
						  sapply(1:nrow(beta), function(x) median(beta[x,(1:numsets-1)*3+2])),
						  sapply(1:nrow(beta), function(x) median(beta[x,(1:numsets-1)*3+3])))
	} else if('betaPH' %in% names(pdata.summary)){
		beta.full = rbind(cbind(sapply(1:nrow(beta), function(x) median(beta[x,(1:numsets-1)*3+1])),
								sapply(1:nrow(beta), function(x) median(beta[x,(1:numsets-1)*3+2])),
								sapply(1:nrow(beta), function(x) median(beta[x,(1:numsets-1)*3+3]))),
						  cbind(sapply(1:nrow(betanph), function(x) median(betanph[x,(1:numsets-1)*3+1])),
								sapply(1:nrow(betanph), function(x) median(betanph[x,(1:numsets-1)*3+2])),
								sapply(1:nrow(betanph), function(x) median(betanph[x,(1:numsets-1)*3+3]))))
	}
	if('k' %in% names(pdata.summary)){
		k.full = cbind(sapply(1:nrow(k), function(x) median(k[x, (1:numsets-1)*3+1])),
					   sapply(1:nrow(k), function(x) median(k[x, (1:numsets-1)*3+2])),
					   sapply(1:nrow(k), function(x) median(k[x, (1:numsets-1)*3+3])))
	}
	if('gamma' %in% names(pdata.summary)){
		gamma.full = cbind(sapply(1:nrow(gamma), function(x) median(gamma[x, (1:numsets-1)*3+1])),
						   sapply(1:nrow(gamma), function(x) median(gamma[x, (1:numsets-1)*3+2])),
						   sapply(1:nrow(gamma), function(x) median(gamma[x, (1:numsets-1)*3+3])))
	}
	
	############ Create the summary table #################
	colnames(hazrate.full) = colnames(pdata.summary$hazardRate)
	rownames(hazrate.full) = rownames(pdata.summary$hazardRate)
	colnames(survrate.full) = colnames(pdata.summary$SurvivalCurve)
	rownames(survrate.full) = rownames(pdata.summary$SurvivalCurve)
	colnames(cumulhaz.full) = colnames(pdata.summary$CumulativeHazard)
	rownames(cumulhaz.full) = rownames(pdata.summary$CumulativeHazard)
	colnames(d.full) = colnames(pdata.summary$d)
	rownames(d.full) = rownames(pdata.summary$d)
	colnames(H.full) = colnames(pdata.summary$H)
	rownames(H.full) = rownames(pdata.summary$H)
	colnames(Rmp.full) = colnames(pdata.summary$Rmp)
	rownames(Rmp.full) = rownames(pdata.summary$Rmp)
	if('beta' %in% names(pdata.summary)){
		colnames(beta.full) = colnames(pdata.summary$beta)
		rownames(beta.full) = rownames(pdata.summary$beta)
	} else if('betaPH' %in% names(pdata.summary)){
		colnames(beta.full) = colnames(pdata.summary$betaPH)
		rownames(beta.full) = c(rownames(pdata.summary$betaPH), rownames(pdata.summary$betaNPH))
	}

	if(GR == TRUE){
		if(length(unique(burnIn)) > 1){
			warning("Burn-in number is not equal across all data sets.  Gelman-Rubin diagnostic test
				 cannot be performed")
			gr = NULL
		}
		if(length(unique(numIters)) > 1){
			warning("Number of MCMC iterations is not equal across all data sets.  Gelman-Rubin diagnostic test
					cannot be performed")
			gr = NULL
		} else {
			gr = gelman.scale(gr.datalist, number.inchain = (numIters[1]-burnIn[1])/thinval, numsets = numsets)
		}
	} 
	output = list(hazardRate = hazrate.full)
	if('beta' %in% names(pdata.summary)){	output = c(output, list(beta = beta.full))	}
	output = c(output, list(SurvivalCurve = survrate.full, CumulativeHazard = cumulhaz.full, d = d.full, 
							H = H.full, Rmp = Rmp.full))
	if('gamma' %in% names(pdata.summary)){	output = c(output, list(gamma = gamma.full))	}
	if('k' %in% names(pdata.summary)){	output = c(output, list(k = k.full))	}
	if(GR == TRUE){	output = c(output, list(gelman.rubin = gr))	}
	
	return(output)
}

##########################################################################################
# The DIC function calculates DIC using the log likelihood.
##########################################################################################
DIC = function(mrhobject, n){
	
	if(is.MRH(mrhobject) == FALSE){ stop("Object is not an MRH object")	}
	if(is.null(colnames(mrhobject))){ if(names(mrhobject)[1] == 'summary'){
		mrhobject = MRH(read.table(paste(mrhobject$outfolder, '/MCMCchains.txt', sep = ''), header = TRUE))
	}}
	nonph.model = TRUE
	if(substr(colnames(mrhobject)[2], 3, 3) == ""){	nonph.model = FALSE	}
	mrhobject = mrhobject[,-1]
	
	if(nonph.model == FALSE){		
		M = log(length(which(substr(colnames(mrhobject), 1, 1) == 'd')))/log(2)
		numhazards = 1
	} else {
		firstgroup = substr(colnames(mrhobject)[2], 4, 4)
		dnames = colnames(mrhobject)[which(substr(colnames(mrhobject), 1, 1) == 'd')]
		numbins = length(which(unlist(strsplit(dnames, '_'))[1:length(dnames)*2] == firstgroup))
		M = log(numbins)/log(2)
		hazgroups = unique(unlist(strsplit(dnames, '_'))[1:length(dnames)*2])
		numhazards = length(hazgroups)
	}
	
	numparams = ncol(mrhobject) - 1 - 2^M*numhazards - 2^M*(numhazards-1)
	Rmp.index = which(substr(colnames(mrhobject), 1, 3) == 'Rmp')
	pruned.Rmps = which(colMeans(mrhobject[,Rmp.index]) == 0.5 &
						sapply(Rmp.index, function(x) var(mrhobject[,x])) == 0)
	numparams = numparams - length(pruned.Rmps)

	neg2loglik = mrhobject[,ncol(mrhobject)]
	neg2loglik.summ = as.data.frame(matrix(summary(neg2loglik), ncol = 1))
    row.names(neg2loglik.summ) = names(summary(neg2loglik))
    names(neg2loglik.summ) = 'value'
    
	DIC = .5*var(neg2loglik) + mean(neg2loglik)
	AIC = max(2*numparams + neg2loglik)
	BIC = NA
	if(!missing(n)){
		BIC = max(neg2loglik + numparams*log(n))
	} else {
		warning("Need number of subjects in data set (n) to calculate BIC. This number can be found in the MCMCInfo.txt file in the output folder.")
	}
	ICtable = as.data.frame(matrix(c(DIC, AIC, BIC), ncol = 1))
	row.names(ICtable) = c('DIC', 'AIC', 'BIC')
    names(ICtable) = 'value'
	return(list(neg2loglik.summ = neg2loglik.summ, ICtable = ICtable))
	
}		

##########################################################################################
# The CalcFunction function calculates the hazard rate, survival rate, and/or 
#	cumulative hazard rate given an MRH object.
##########################################################################################
CalcFunction = function(mrhobject, function.type = c('h', 'H', 'S'), maxStudyTime, alpha.level = 0.05){

	censortime = NULL
	if(!missing(maxStudyTime)){	censortime = maxStudyTime	}
	################ Find M, get the chains #############
	if(is.MRH(mrhobject) == FALSE){ stop("Object is not an MRH object")	}
	if(is.null(colnames(mrhobject))){ if(names(mrhobject)[1] == 'summary'){
		censortime = mrhobject$maxStudyTime
		mrhobject = MRH(read.table(paste(mrhobject$outfolder, '/MCMCchains.txt', sep = ''), header = TRUE))
	}}
	if(is.null(censortime) & 'h' %in% function.type){ stop("Maximum study time needed for hazard rate calculation")	}
	nonph.model = TRUE
	if(substr(colnames(mrhobject)[2], 3, 3) == ""){	nonph.model = FALSE	}
	mrhobject = mrhobject[,-1]
	
	if(nonph.model == FALSE){		
		M = log(length(which(substr(colnames(mrhobject), 1, 1) == 'd')))/log(2)
		numhazards = 1
	} else {
		firstgroup = substr(colnames(mrhobject)[2], 4, 4)
		dnames = colnames(mrhobject)[which(substr(colnames(mrhobject), 1, 1) == 'd')]
		numbins = length(which(unlist(strsplit(dnames, '_'))[1:length(dnames)*2] == firstgroup))
		M = log(numbins)/log(2)
		hazgroups = unique(unlist(strsplit(dnames, '_'))[1:length(dnames)*2])
		numhazards = length(hazgroups)
	}
	q.probs = c(alpha.level/2, .5, 1-alpha.level/2)
	if('h' %in% function.type){	binwidth = censortime/2^M	}
	
	############# Calculate hazard rate ##############
	if('h' %in% function.type){
		hazrate = sapply(1:(2^M*numhazards), function(x) quantile(mrhobject[,x], probs = q.probs))/binwidth
		dnames = paste('h.bin', 1:2^M, sep = '')
		if(numhazards > 1){
			dnames = paste(rep(dnames, numhazards), rep(paste('.group', 1:numhazards, sep = ''), each = 2^M), sep = '')
		}
		colnames(hazrate) = dnames
		rownames(hazrate) = paste('q', c(alpha.level/2, .5, 1-alpha.level/2), sep = '.')
		output = list(hazrate = hazrate)
	} 
	############# Calculate the cumulative hazard and/or survival functions ##############
	if('H' %in% function.type || 'S' %in% function.type){
		mcmc.cumulhaz = NULL
		for(hazctr in 1:numhazards){
			mcmc.cumulhaz = cbind(mcmc.cumulhaz, matrix(0, ncol = 1, nrow = nrow(mrhobject)))
			for(colctr in 1:2^M){
				mcmc.cumulhaz = cbind(mcmc.cumulhaz, mcmc.cumulhaz[,colctr+(hazctr-1)*(2^M+1)]+mrhobject[,colctr+(hazctr-1)*2^M])
			}
		}
		if('H' %in% function.type){
			cumulhaz = sapply(1:((2^M+1)*numhazards), function(x) quantile(mcmc.cumulhaz[,x], probs = q.probs))
			cumulnames = c('H.start', paste('H.bin', 1:2^M, sep = ''))
			if(numhazards > 1){
				cumulnames = paste(rep(cumulnames, numhazards), rep(paste('.group', 1:numhazards, sep = ''), each = 2^M+1), sep = '')
			}
			colnames(cumulhaz) = cumulnames
			rownames(cumulhaz) = paste('q', c(alpha.level/2, .5, 1-alpha.level/2), sep = '.')
			if('h' %in% function.type){	output = c(output, list(cumulhaz = cumulhaz))
			} else {	output = list(cumulhaz = cumulhaz)	}
		}
		if('S' %in% function.type){
			mcmc.St = exp(-mcmc.cumulhaz)
			survfxn = sapply(1:((2^M+1)*numhazards), function(x) quantile(mcmc.St[,x], probs = q.probs))
			survnames = c('S.start', paste('S.bin', 1:2^M, sep = ''))
			if(numhazards > 1){
				survnames = paste(rep(survnames, numhazards), rep(paste('.group', 1:numhazards, sep = ''), each = 2^M+1), sep = '')
			}
			colnames(survfxn) = survnames
			rownames(survfxn) = paste('q', c(alpha.level/2, .5, 1-alpha.level/2), sep = '.')
			if(length(function.type) > 1){ output = c(output, list(survfunction = survfxn))
			} else {	output = list(survfunction = survfxn)	}
		}
	}
	return(output)
}







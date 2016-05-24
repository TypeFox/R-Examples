###########################################################################################
# Contains the functions needed to: convert to an MRH object, summary, and plot for PH model.
###########################################################################################
library(survival)
library(coda)

MRH = function(x) UseMethod("MRH")

MRH.default = function(x){
	if(nrow(x) == 1 | ncol(x) < 5){ 
		stop("Data set provided is not an MRH MCMC chain")
	}
	x = as.matrix(x)
	class(x) <- "MRH"
	x
}

is.MRH = function(x){	if(inherits(x, "MRH")){	TRUE } else {	FALSE	}	}
as.MRH = function(x){	UseMethod("MRH")	}

# print will eventually contain the call, etc (see lm results)
summary.MRH = function(object, alpha.level = .05, maxStudyTime,...){
	
	if(is.null(colnames(object))){ if(names(object)[1] == 'summary'){
		maxStudyTime = object$maxStudyTime
		object = MRH(read.table(paste(object$outfolder, '/MCMCchains.txt', sep = ''), header = TRUE))
	}} 
	if(missing(maxStudyTime)){ stop("Maximum study time (maxStudyTime) needed for hazard rate calculation. The maximum study time can be found in the MCMCInfo.txt file in the output folder.")	}
	
	substr.5 = substr(colnames(object), 1, 5)
	substr.4 = substr(colnames(object), 1, 4)
	substr.3 = substr(colnames(object), 1, 3)
	substr.1 = substr(colnames(object), 1, 1)
	
	any.betas = any.ks = any.gammas = FALSE
	if(length(which(substr.5 == 'gamma')) > 0){
		any.gammas = TRUE
		gamma.loc = which(substr.5 == 'gamma')
	}
	if(length(which(substr.4 == 'beta')) > 0){	
		any.betas = TRUE	
		beta.loc = which(substr.4 == 'beta')
	}
	if(length(which(substr.1 == 'k')) > 0){	
		any.ks = TRUE	
		k.loc = which(substr.1 == 'k')
	}
	H.loc = which(substr.1 == 'H')
	Rmp.loc = which(substr.3 == 'Rmp')
	d.loc = which(substr.1 == 'd')
	
	title.lb = alpha.level/2*1000
	if(alpha.level/2*1000 < 100){	title.lb = paste('0', title.lb, sep = '')	}
	if(alpha.level/2*1000 < 10){	title.lb = paste('0', title.lb, sep = '')	}
	title.ub = (1-alpha.level/2)*1000
	
	################ Calculate the outputs for each parameter ################
	### Calculate hazard rate, survival function, and cumulative hazard
	allfunctions = CalcFunction(as.MRH(object), maxStudyTime = maxStudyTime, alpha.level = alpha.level)
	hroutput = as.data.frame(t(allfunctions$hazrate)[,c(2,1,3)])
	names(hroutput) = c('hrEst', paste('hrq', title.lb, sep = '.'), paste('hrq', title.ub, sep = '.'))
	row.names(hroutput) = colnames(allfunctions$hazrate)
	
	survoutput = as.data.frame(t(allfunctions$survfunction)[,c(2,1,3)])
	names(survoutput) = c('SEst', paste('Sq', title.lb, sep = '.'), paste('Sq', title.ub, sep = '.'))
	row.names(survoutput) = colnames(allfunctions$survfunction)
	
	cumulhoutput = as.data.frame(t(allfunctions$cumulhaz)[,c(2,1,3)])
	names(cumulhoutput) = c('HEst', paste('Hq', title.lb, sep = '.'), paste('Hq', title.ub, sep = '.'))
	row.names(cumulhoutput) = colnames(allfunctions$cumulhaz)
	
	### Calculate beta if there are any ###
	if(any.betas == TRUE){
		betas = t(sapply(beta.loc, function(x) 
						 quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2))))
		betaoutput = as.data.frame(betas)
		names(betaoutput) = c('betaEst', paste('betaq', title.lb, sep = '.'), paste('betaq', title.ub, sep = '.'))
		row.names(betaoutput) = colnames(object)[beta.loc]
	}
	
	### Calculate d
	dInts = t(sapply(d.loc, function(x) quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2))))
	doutput = as.data.frame(dInts)
	names(doutput) = c('dEst', paste('dq', title.lb, sep = '.'), paste('dq', title.ub, sep = '.'))
	row.names(doutput) = colnames(object)[d.loc]
	
	#### Calculate the Rmps 
	RmpInts = as.data.frame(t(sapply(Rmp.loc, 
									 function(x) quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2)))))
	names(RmpInts) = c('RmpEst', paste('Rmpq', title.lb, sep = '.'), paste('Rmpq', title.ub, sep = '.'))
	row.names(RmpInts) = colnames(object)[Rmp.loc]
	
	#### Calculate the cumulative hazard H
	Houtput = as.data.frame(t(sapply(H.loc, function(x) quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2)))))
	names(Houtput) = c('HEst', paste('Hq', title.lb, sep = '.'), paste('Hq', title.ub, sep = '.'))
	row.names(Houtput) = colnames(object)[H.loc]
	
	#### Calculate gamma if it exists ###
	if(any.gammas == TRUE){
		gammaoutput = as.data.frame(t(sapply(gamma.loc, function(x)
                                    quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2)))))
		names(gammaoutput) = c('gammampEst', paste('gammampq', title.lb, sep = '.'), paste('gammampq', title.ub, sep = '.'))
		row.names(gammaoutput) = colnames(object)[gamma.loc]
	}
	
	#### Calculate k if it exists ###
	if(any.ks == TRUE){
		koutput = as.data.frame(t(sapply(k.loc,
                                function(x) quantile(object[,x], probs = c(.5, alpha.level/2, 1-alpha.level/2)))))
		names(koutput) = c('kEst', paste('kq', title.lb, sep = '.'), paste('kq', title.ub, sep = '.'))
		row.names(koutput) = colnames(object)[k.loc]
	}
	
	output = list(hazardRate = hroutput)
	if(any.betas == TRUE){	output = c(output, list(beta = betaoutput))	}
	output = c(output, list(SurvivalCurve = survoutput, CumulativeHazard = cumulhoutput, 
							d = doutput, H = Houtput, Rmp = RmpInts))
	
	if(any.gammas == TRUE){	output = c(output, list(gamma = gammaoutput))	}
	if(any.ks == TRUE){	output = c(output, list(k = koutput))	}
	
	return(output)
}

plot.MRH = function(x, maxStudyTime, main = "", xlab = 'Time', ylab = 'Hazard Rate', plot.type = 'h', interval = TRUE,
alpha.level = 0.05, smooth.graph = FALSE, smooth.df = NULL, combine.graphs = TRUE, log.ratio = TRUE, ...){
	
	censortime = NULL
	if(!missing(maxStudyTime)){	censortime = maxStudyTime	}
	if(is.MRH(x) == FALSE){ stop("Object is not an MRH object")	}
	if(is.null(colnames(x))){ if(names(x)[1] == 'summary'){
		censortime = x$maxStudyTime
		x = MRH(read.table(paste(x$outfolder, '/MCMCchains.txt', sep = ''), header = TRUE))
	}} 
	if(is.null(censortime) & 'h' == plot.type){ stop("Maximum study time needed for hazard rate calculation")	}

	if(substr(colnames(x)[2], 3, 3) == ""){
		callPHplot(x, censortime, main = main, xlab = xlab, ylab = ylab, plot.type = plot.type, 
				   smooth.graph = smooth.graph, smooth.df = smooth.df, conf.int = interval, alpha.level = alpha.level)
	} else {
		callBAplot(x, censortime, main = main, xlab = xlab, ylab = ylab, plot.type = plot.type, 
				   smooth.graph = smooth.graph, smooth.df = smooth.df, combine.graphs = combine.graphs, 
				   conf.int = interval, alpha.level = alpha.level, log.ratio = log.ratio)
	}
}

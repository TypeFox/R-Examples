######################################################################################
# Created 23Jan12: This file contains the functions that are needed to determine
#	"convergence".  Includes:
#	- thinAmt: Tests the autocorrelation to make sure there is enough thinning.  
#				Based on code found here: http://stats.stackexchange.com/questions/4258/can-i-semi-automate-mcmc-convergence-diagnostics-to-set-the-burn-in-length
#	- convergeFxn: Runs the geweke.diag test and counts how many parameters are over
#				z(alpha = .01/2)
#	- convergenceGraphs: XXX
######################################################################################

thinAmt = function(thinorig, dataset){

	lagPoss = c(0, 1, 5, 10, 15)
	acm = autocorr.diag(dataset, lags = lagPoss)
	
	# This line of code find the first place/lag where autocorrelation disappears
	# (is calculated negative) for each variable, and then the minimum of 
	# those lags is the new thinning amount.
	# - acm < 0 makes a matrix of TRUE/FALSE for values less than 0
	#	'2' in the apply function tells apply to go by column
	# - match() finds the first TRUE in each column, which denotes 
	#	the first place autocorrelation disappears?
	minThin = min(apply(acm < 0.1, 2, function(x) match(TRUE, x)), na.rm = TRUE)
	
	returnthin = lagPoss[minThin]*thinorig
	
	return(returnthin)
}
	
convergeFxn = function(dataset){

	# Geweke diagnostic first #
	zvals = geweke.diag(dataset)$z
	if(length(which(abs(zvals) > qnorm(.99))) > 0){	
		return(FALSE)
	} else {	
		# Heidelberg-Welch if Geweke passes #
		HWvals = heidel.diag(dataset)[1:ncol(dataset),]
		if(length(which(HWvals[,1] == 0)) > 0){
		return(FALSE)
		} else {
			return(TRUE)	
		}
	}
}

convergenceGraphs = function(dataset, graphname, burn.in, thin.val){

    paramNames = names(dataset)
	byThrees = seq(1, ncol(dataset)-3, by = 3)

	masize = floor(nrow(dataset)/10)
	if(masize == 0){	masize = 1	}
	
	for(graphCtr in byThrees){
		
		pdf(file = paste(graphname, (graphCtr-1)/3+1, '.pdf', sep = ''), width = 8, height = 5)
		par(mfrow = c(3,4))

		subGraphFxn(dataset[,graphCtr], paramNames[graphCtr], burn.in, thin.val, masize)  
		subGraphFxn(dataset[,graphCtr+1], paramNames[graphCtr+1], burn.in, thin.val, masize)
		subGraphFxn(dataset[,graphCtr+2], paramNames[graphCtr+2], burn.in, thin.val, masize)
		dev.off()
	}
	pdf(file = paste(graphname, (graphCtr+2)/3+1, '.pdf', sep = ''), width = 8, height = 5)
	par(mfrow = c(3,4))
	for(graphCtr in (max(byThrees)+3):ncol(dataset)){
		subGraphFxn(dataset[,graphCtr], paramNames[graphCtr], burn.in, thin.val, masize)
	}
	dev.off()
}
			
subGraphFxn = function(vectorVals, titleName, burninval, thinval, MAsize){

	plot((1:length(vectorVals))*thinval+burninval, vectorVals, 
		 main = paste('Trace', titleName), xlab = 'Iteration', 
		 ylab = 'Parameter value', type = 'l')

	plot(density(vectorVals), main = paste('Density', titleName), xlab = '')

	plot((1:(length(vectorVals)-(MAsize-1)))*thinval+burninval, MA(vectorVals, MAsize), 
		 type = 'l', main = paste('MA', titleName), 
		 xlab = 'Iteration', ylab = 'Average')

	acvals = autocorr(mcmc(vectorVals, thin = thinval, start = burninval), lags = 0:20)[,,1]
	plotlags = (0:20)*thinval
	plot(rep(plotlags[1], 2), c(0, acvals[1]), type = 'l', ylim = c(-1,1), xlim = range(plotlags),
		 main = paste('Autocorrelation', titleName), xlab = 'Lag', ylab = 'AC')

	for(i in 2:length(plotlags)){
		points(rep(plotlags[i], 2), c(0, acvals[i]), type = 'l')
	}
}

MA = function(vectorVals, avgSize, numsets){
	
	sapply(1:(length(vectorVals)-(avgSize-1)), function(x) sum(vectorVals[x:(x+(avgSize-1))]))/avgSize
}

gelman.scale = function(datalist, number.inchain, numsets){
	
	# Find the variances for each parameter in each chain.
	sj_sq = meanj = NULL
	for(j in 1:numsets){
		sj_sq = rbind(sj_sq, sapply(1:ncol(datalist[[j]]), function(x) var(datalist[[j]][,x])))
		meanj = rbind(meanj, colMeans(datalist[[j]]))
	}
	
	# W is the mean of the variances for each parameter	
	W = colMeans(sj_sq)
	# B is the weighted squared difference between means for each parameter
	B = number.inchain*sapply(1:ncol(meanj), function(x) var(meanj[,x]))
	
	# Var theta is the variance of the stationary distribution as a wtd avg of W and B
	VarThetas = (1-1/number.inchain)*W + 1/number.inchain*B
	
	# R is the scale reduction factor
	R = sqrt(VarThetas/W)
	
	output = as.data.frame(round(R, 2))
	row.names(output) = colnames(datalist[[1]])
	colnames(output) = 'Scale Reduction Factor'
	
	return(output)
}
	


######################################################################################
# Created 18Jan12: This file contains an assortment of functions that are needed to 
#	run the Gibbs sampler, but the user will not have use for them.  
# This file contains:
# calc_dfast - a quick way to calculate the ds.
# calc_mat01 - calculates the 0/1 matrix used in calculating d.
# getbin.and.ratio - For each person, this gets the bin they failed in and the 
#	ratio of time they spent in that bin.
# calcXIndics - Used in the computation of the Rmp posteriors.
# AnalyzeComplete - Analyzes the final parameter estimates from the Gibbs sampler by
#	calculating estimates and credible intervals for d, beta, H, and Rmp, and plotting
#	the results.
# PlotCIs - Plots the estimates and credible intervals (used by AnalyzeComplete)
# PlotCIsSmooth - Plots the smoothed estimates and credible intervals (used by AnalyzeComplete)
######################################################################################
library(survival)
library(coda)
###########################################################################################
# The calc_dfast function calculates d faster given an M, Rmp values,  H00, and the 
#	0/1 matrix is passed to it.
###########################################################################################
calc_dfast = function(M, Rmpvals, H00val, mat01){
	
	Rmpvec = rep(NA, length = sum(2^(1:M)))
	Rmpvec[1:(sum(2^(1:M))/2)*2-1] = Rmpvals
	Rmpvec[1:(sum(2^(1:M))/2)*2] = 1-Rmpvals
	ds = H00val*exp(mat01%*%log(Rmpvec))
	
	return(ds)
}
###########################################################################################
# The calc_mat01 function is used to calculate the matrix needed to calculate the ds 
#	and other parameters (such as F).  The user must enter M, the highest number of 
#	levels for the tree.
###########################################################################################
calc_mat01 = function(M){
	
	mat_01 = matrix(0, ncol = sum(2^(1:M)), nrow = 2^M)
	
	for(i in 1:M){
		incr.unit = 2^(M-i)
		
		for(j in 0:(2^i-1)){
			mat_01[(j*incr.unit+1):((j+1)*incr.unit), sum(2^(0:(i-1)))+j] = 1
		}
	}
	
	return(mat_01)
}

###########################################################################################
# The get.bin.and.ratio function calculates the correct bin and ratio for each Ti.  This 
#	is done by first calculating the interval length, which is equal to the censor time (TJ)
#	divided by the total number of bins 2^M.  The bin number is the number of bins that 
#	Ti goes past, and then the ratio is the length of the interval the Ti falls into,
#	divided by the interval length.  For example, if Ti = .75, 1.5, 2.3, and 4.2, and the 
#	interval length is 1, then the bin number would be 0, 1, 2, 4, and the ratio would be
#	.75, .50, .30, .20.  The user must enter the censoring time (TJ), the observed 
#	failure times for each subject (Ti), and the maximum hierarchy level (M).
###########################################################################################
getbin.and.ratio = function(TJ, Ti, M){
	
	bin = rep(0, length(Ti))
	interval.length = TJ/2^M
	ratio = rep(NA, length(Ti))
	
	for(i in 1:(2^M)){
	# For subjects who have a Ti greater than or equal to the current time (time is 
	# calculated as i*interval.length), put them in bin i.	
		bin[which(Ti >= i*interval.length)] = i
		# For subjects who have a Ti less than the current time, but greater than or equal
		# to the previous time, calculate the ratio.
		ratio[which((i-1)*interval.length <= Ti & Ti < i*interval.length)] = 
		(Ti[which((i-1)*interval.length <= Ti & Ti < i*interval.length)]-(i-1)*interval.length)/interval.length
	}
	ratio[which(Ti >= 2^M*interval.length)] = 1
	
	output = as.data.frame(cbind(bin, ratio))
	names(output) = c('bin','ratio')
	return(output)
}

###########################################################################################
# This function is used to calculate the necessary h(Ti) values used in the posterior of
#	Rmp.  Because h(Ti) is the dj that the Ti falls into, and because the djs are 
#	functions of the Rmps, the Rmps that do not actually belong to the selected m and p
#	need to be eliminated.  What this reduces to is for each dj and Rmp (for a given m and
#	p), the dj either belongs to Rmp, 1-Rmp or neither.  This code denotes, for each person,
#	whether the person will contribute an Rmp, a 1-Rmp, or neither to the likelihood fxn.
# NOTE: This only needs to be called ONCE for estimating the parameters of a data set. 
###########################################################################################
calcXIndics = function(mat01, bin_ratio, M){
	
	# xIndic will hold the indicators for Rmp for each for each person depending 
	# on which dj the Ti falls into.
	# xIndic will hold the indicators for 1-Rmp for each for each person depending 
	# on which dj the Ti falls into.
	xIndic = one_xIndic = matrix(NA, ncol = sum(2^(1:M))/2, nrow = nrow(bin_ratio))
	di.num = bin_ratio$bin+1
	di.num[which(di.num > 2^M)] = 2^M
	for(m in 1:(2^M)){
		if(length(which(di.num == m)) > 0){
			# For subjects with Ti in bin m, give the dj indicators with j = m #
			xIndic[which(di.num == m),] = matrix(c(mat01[m, 1:(sum(2^(1:M))/2)*2-1]), 
							nrow = length(which(di.num == m)), byrow = T, ncol = 2^M-1)
			one_xIndic[which(di.num == m),] = matrix(c(mat01[m, 1:(sum(2^(1:M))/2)*2]),
							nrow = length(which(di.num == m)), byrow = T, ncol = 2^M-1)
		}
	}
	
	return(cbind(xIndic, one_xIndic))
}

###########################################################################################
# The AnalyzeComplete function does the following: 
#	1. Calculate the credible intervals for d, beta, H, and Rmp.  
#	2. Plot the hazard rate with the credible intervals.
#	3. Plot the smoothed hazard rate with credible intervals.
# The graphs are saved to files, and the estimates and intervals are returned to the user.
###########################################################################################
AnalyzeComplete = function(M, censortime, finaldataset, outfilename, graphs = TRUE, numhazards, namesHazGroups){
	
	# Get the summaries for the data set
	param.summs = summary(as.MRH(finaldataset), maxStudyTime = censortime)

    ############# Calculate and graph the hazard rate ############
	if(graphs == TRUE){
		for(hazCtr in 1:numhazards){
			# Plot the ds over time
			if(numhazards == 1){
				plotCIs(M, censortime, param.summs$hazardRate[1:2^M+(hazCtr-1)*2^M,c(2,1,3)], 
						paste(outfilename, 'HazardRateGraph', sep = ''))
			} else {
				plotCIs(M, censortime, param.summs$hazardRate[1:2^M+(hazCtr-1)*2^M,c(2,1,3)], 
						paste(outfilename, 'HazardRateGraph', namesHazGroups[hazCtr], sep = ''))
			}
			if(M > 1){
				# Calculate and plot the smooth ds over time 
				timePts = seq(censortime/2^M, censortime, by = censortime/2^M)-.5*censortime/2^M
				smoothdInts = smooth.spline(param.summs$hazardRate[1:2^M+(hazCtr-1)*2^M,1]~timePts, df = 2^M/2)$y
				smoothdInts = cbind(smoothdInts, smooth.spline(param.summs$hazardRate[1:2^M+(hazCtr-1)*2^M,2]~timePts, df = 2^M/2)$y)
				smoothdInts = cbind(smoothdInts, smooth.spline(param.summs$hazardRate[1:2^M+(hazCtr-1)*2^M,3]~timePts, df = 2^M/2)$y)
				if(numhazards == 1){
					plotCIsSmooth(M, censortime, smoothdInts, paste(outfilename, 'HazardRateSmoothGraph', sep = ''))
				} else {
					plotCIsSmooth(M, censortime, smoothdInts, paste(outfilename, 'HazardRateSmoothGraph', namesHazGroups[hazCtr], sep = ''))
				}
			}
		}
	}
	
	############ Create the summary table #################
	summaryTable = ACniceOutput(param.summs$hazardRate)
    if('beta' %in% names(param.summs)){    any.betas = TRUE    } else {    any.betas = FALSE   }
	if(any.betas == TRUE){ summaryTable = rbind(summaryTable, ACniceOutput(param.summs$beta))	}
		
	# Organize the summary table
	summaryTable = as.data.frame(summaryTable)
	names(summaryTable) = c('Estimate', 'CredibleInterval')
	temprownames =  rownames(param.summs$hazardRate)
	if(any.betas == TRUE){	temprownames = c(temprownames, rownames(param.summs$beta))	}
    row.names(summaryTable) = temprownames
	
	output = c(list(summary = summaryTable), param.summs)
	
	return(output)
}

ACniceOutput = function(output){
	
	output2 = cbind(round(output[,1], 3), paste("(",round(output[,2], 3), ", ", round(output[,3], 3), ")", sep = ''))
		return(output2)
}
																							
GraphNPbetas = function(M, dests, np.betaests, numhazgroups, hazgroupnames, censortime, file){

	##################### RAW DS, SEPARATE GRAPHS ##############
	xvals = pbfxn(2^M, censortime/2^M, dests[1:2^M,1])$x
	binwidth = censortime/2^M
	pdf(file = paste(file, 'SeparateHazardRates.pdf', sep  = ''), width = 7, height = 2.5*ceiling(numhazgroups/2))
	par(mfrow = c(ceiling(numhazgroups/2), 2), mai = c(.35, .1, .25, 0), oma = c(0, 2, 2, 1), tck = -.02)
	for(hazCtr in 1:numhazgroups){
		plot(xvals, pbfxn(2^M, censortime/2^M, dests[1:2^M+2^M*(hazCtr-1),1]/binwidth)$y, type = 'l', lwd = 3, 
			 ylab = '', xlab = '', ylim = range(dests/binwidth), 
			 main = paste(hazgroupnames[hazCtr], 'group'), axes = FALSE)
		points(xvals, pbfxn(2^M, censortime/2^M, dests[1:2^M+2^M*(hazCtr-1),2]/binwidth)$y, type = 'l', lty = 2)
		points(xvals, pbfxn(2^M, censortime/2^M, dests[1:2^M+2^M*(hazCtr-1),3]/binwidth)$y, type = 'l', lty = 2)
		box()
		axis(1, at = round(xvals), labels = round(xvals), cex.axis = .75, padj = -2)
		mtext('Time', side = 1, padj = 1.5, cex = .8)
		if(hazCtr %% 2 != 0){	axis(2, padj = 1.25)	}
	}
	mtext('Hazard Rates by Group', outer = TRUE, cex = 1.4)
	dev.off()
	
	##################### RAW DS, COMBINED ##############
	pdf(file = paste(file, 'HazardRates.pdf', sep  = ''))
	plot(xvals, pbfxn(2^M, censortime/2^M, dests[1:2^M,1]/binwidth)$y, type = 'l', lwd = 2,
		 ylab = 'Time', xlab = 'Hazard rate', main = 'Hazard Rates', ylim = range(dests[,1]/binwidth))
	for(hazCtr in 2:numhazgroups){
		points(xvals, pbfxn(2^M, censortime/2^M, dests[1:2^M+2^M*(hazCtr-1),1]/binwidth)$y, 
			   type = 'l', lwd = 2, col = hazCtr)
	}
	legend(x = "topright", paste('group', hazgroupnames), fill = 1:numhazgroups, cex = 1.2, 
		   inset = max(c(min(dests[,1]/binwidth), .01))/5)
	dev.off()
	
	##################### RAW BETAS, SEPARATE GRAPHS ##############
	widthval = 7
	if(numhazgroups > 2){	
		plotnums = c(ceiling((numhazgroups-1)/2), 2)	
		heightval = 3.5*ceiling((numhazgroups-1)/2)
	} else {	
		plotnums = c(1,1)	
		heightval = widthval
	}
	pdf(file = paste(file, 'SeparateHazardRatios.pdf', sep = ''), height = heightval, width = widthval)
	par(mfrow = plotnums, mai = c(.35, .1, .25, 0), oma = c(0, 2, 2, 1), tck = -.02)
	for(ratioCtr in 1:(numhazgroups-1)){
		plot(c(0, censortime), c(0, 0), type = 'l', col = 'grey', lwd = 2, ylab = '', 
			 ylim = range(np.betaests), main = paste('Group', hazgroupnames[ratioCtr+1]), axes = FALSE)
		points(xvals, pbfxn(2^M, censortime/2^M, np.betaests[1:2^M+2^M*(ratioCtr-1),1])$y, type = 'l', lwd = 3)
		points(xvals, pbfxn(2^M, censortime/2^M, np.betaests[1:2^M+2^M*(ratioCtr-1),2])$y, type = 'l', lty = 2)
		points(xvals, pbfxn(2^M, censortime/2^M, np.betaests[1:2^M+2^M*(ratioCtr-1),3])$y, type = 'l', lty = 2)
		box()
		axis(1, at = round(xvals), labels = round(xvals), cex.axis = .75, padj = -2)
		mtext('Time', side = 1, padj = 1.5, cex = .8)
		if(ratioCtr %% 2 != 0){	axis(2, padj = 1.25)	}
	}
	mtext(paste('Hazard Ratios: Comparison to Group', hazgroupnames[1]), outer = TRUE, cex = 1.4)
	dev.off()
	
	##################### RAW BETAS, ONE GRAPH ##############
	if(numhazgroups > 2){
		pdf(file = paste(file, 'HazardRatios.pdf', sep  = ''))
		plot(c(0, censortime), c(0, 0), type = 'l', col = 'grey', lwd = 2, 
			 ylab = 'Time', xlab = 'Hazard ratio', 
			 main = paste('Hazard Ratios: Comparison to Group', hazgroupnames[1]), ylim = range(np.betaests[,1]))
		for(ratioCtr in 2:numhazgroups-1){
			points(xvals, pbfxn(2^M, censortime/2^M, np.betaests[1:2^M+2^M*(ratioCtr-1),1])$y, 
				   type = 'l', lwd = 2, col = ratioCtr)
		}
		legend(x = "topleft", paste('group', hazgroupnames[2:numhazgroups]), fill = 2:numhazgroups-1, cex = 1.2, 
			   inset = max(c(min(np.betaests[,1]), .01))/5)
		dev.off()
	}
	
	##################### SMOOTH DS, SEPARATE GRAPHS ##############
	smoothds = smoothbetas = NULL
	timepts = 1:2^M*binwidth-.5*binwidth
	for(hazCtr in 1:numhazgroups){
		smoothds = rbind(smoothds, cbind(smooth.spline(dests[1:2^M+(hazCtr-1)*2^M,1]/binwidth ~ timepts, df = 2^M/2)$y, 
										 smooth.spline(dests[1:2^M+(hazCtr-1)*2^M,2]/binwidth ~ timepts, df = 2^M/2)$y, 
										 smooth.spline(dests[1:2^M+(hazCtr-1)*2^M,3]/binwidth ~ timepts, df = 2^M/2)$y))
	}
	for(hazCtr in 1:(numhazgroups-1)){
		smoothbetas = rbind(smoothbetas, cbind(smooth.spline(np.betaests[1:2^M+(hazCtr-1)*2^M,1] ~ timepts, df = 2^M/2)$y,
											   smooth.spline(np.betaests[1:2^M+(hazCtr-1)*2^M,2] ~ timepts, df = 2^M/2)$y,
											   smooth.spline(np.betaests[1:2^M+(hazCtr-1)*2^M,3] ~ timepts, df = 2^M/2)$y))
	}
	pdf(file = paste(file, 'SeparateHazardRates_smooth.pdf', sep  = ''), 
		width = 7, height = 2.5*ceiling(numhazgroups/2))
	par(mfrow = c(ceiling(numhazgroups/2), 2), mai = c(.35, .1, .25, 0), oma = c(0, 2, 2, 1), tck = -.02)
	for(hazCtr in 1:numhazgroups){
		plot(timepts, smoothds[1:2^M+(hazCtr-1)*2^M,1], type = 'l', lwd = 3, ylab = '', xlab = '', 
			 ylim = range(smoothds), main = paste(hazgroupnames[hazCtr], 'group'), axes = FALSE)
		points(timepts, smoothds[1:2^M+(hazCtr-1)*2^M,2], type = 'l', lty = 2)
		points(timepts, smoothds[1:2^M+(hazCtr-1)*2^M,3], type = 'l', lty = 2)
		box()
		axis(1, at = round(xvals), labels = round(xvals), cex.axis = .75, padj = -2)
		mtext('Time', side = 1, padj = 1.5, cex = .8)
		if(hazCtr %% 2 != 0){	axis(2, padj = 1.25)	}
	}
	mtext('Hazard Rates by Group', outer = TRUE, cex = 1.4)
	dev.off()
	
	##################### SMOOTH DS, ONE GRAPH ##############
	pdf(file = paste(file, 'HazardRates_smooth.pdf', sep  = ''))
	plot(timepts, smoothds[1:2^M,1], type = 'l', lwd = 2, ylab = 'Time', xlab = 'Hazard rate', 
		 main = 'Hazard Rates', ylim = range(smoothds[,1]))
	for(hazCtr in 2:numhazgroups){
		points(timepts, smoothds[1:2^M+(hazCtr-1)*2^M,1], type = 'l', lwd = 2, col = hazCtr)
	}
	legend(x = "topright", paste('group', hazgroupnames), fill = 1:numhazgroups, cex = 1.2, 
		   inset = max(c(min(dests[,1]/binwidth), .01))/5)
	dev.off()

	##################### SMOOTH BETAS, SEPARATE GRAPHS ##############
	widthval = 7
	if(numhazgroups > 2){	
		plotnums = c(ceiling((numhazgroups-1)/2), 2)	
		heightval = 3.5*ceiling((numhazgroups-1)/2)
	} else {	
		plotnums = c(1,1)	
		heightval = widthval
	}
	pdf(file = paste(file, 'SeparateHazardRatios_smooth.pdf', sep = ''), height = heightval, width = widthval)
	par(mfrow = plotnums, mai = c(.35, .1, .25, 0), oma = c(0, 2, 2, 1), tck = -.02)
	for(ratioCtr in 1:(numhazgroups-1)){
		plot(c(0, censortime), c(0, 0), type = 'l', col = 'grey', lwd = 2, ylab = '', 
			 ylim = range(smoothbetas), main = paste('Group', hazgroupnames[ratioCtr+1]), axes = FALSE)
		points(timepts, smoothbetas[1:2^M+2^M*(ratioCtr-1),1], type = 'l', lwd = 3)
		points(timepts, smoothbetas[1:2^M+2^M*(ratioCtr-1),2], type = 'l', lty = 2)
		points(timepts, smoothbetas[1:2^M+2^M*(ratioCtr-1),3], type = 'l', lty = 2)
		box()
		axis(1, at = round(xvals), labels = round(xvals), cex.axis = .75, padj = -2)
		mtext('Time', side = 1, padj = 1.5, cex = .8)
		if(ratioCtr %% 2 != 0){	axis(2, padj = 1.25)	}
	}
	mtext(paste('Hazard Ratios: Comparison to Group', hazgroupnames[1]), outer = TRUE, cex = 1.4)
	dev.off()
	
	##################### SMOOTH BETAS, ONE GRAPH ##############
	if(numhazgroups > 2){
		pdf(file = paste(file, 'HazardRatios_smooth.pdf', sep  = ''))
		plot(c(0, censortime), c(0, 0), type = 'l', col = 'grey', lwd = 2, ylab = 'Time', xlab = 'Hazard ratio', 
			 main = paste('Hazard Ratios: Comparison to Group', hazgroupnames[1]), ylim = range(smoothbetas[,1]))
		for(ratioCtr in 2:numhazgroups-1){
			points(timepts,smoothbetas[1:2^M+2^M*(ratioCtr-1),1], type = 'l', lwd = 2, col = ratioCtr)
		}
		legend(x = "topleft", paste('group', hazgroupnames[2:numhazgroups]), fill = 2:numhazgroups-1, cex = 1.2, 
			   inset = max(c(min(np.betaests[,1]), .01))/5)
		dev.off()
	}
	
}	
	
	
###########################################################################################
# The plotCIs function plots the credible intervals for the actual hazard estimates.  The
#	user must provide M, the censortime, the estimates and their bounds, and the graph name.
#	This function is used in the AnalyzeComplete function.
###########################################################################################
plotCIs = function(M, censortime, plotvalues, graphname, 
mainTitle = 'Hazard rate', xlabel = 'Time', ylabel = 'Hazard estimate'){
	
	xvals = pbfxn(2^M, censortime/2^M, plotvalues[,1])
	pdf(file = paste(graphname, '.pdf', sep = ''), width = 8, height = 6)
	plot(xvals$x, xvals$y, main = mainTitle, xlab = xlabel, ylab = ylabel, type = 'l',
		 lty = 2, ylim = range(plotvalues))
	points(xvals$x, pbfxn(2^M, censortime/2^M, plotvalues[,2])$y, type = 'l', lwd = 3)
	points(xvals$x, pbfxn(2^M, censortime/2^M, plotvalues[,3])$y, type = 'l', lty = 2)
	dev.off()
}

###########################################################################################
# The plotCIsSmooth function plots the credible intervals for the smoothed hazard estimates.  
#	The user must provide M, the censortime, the smoothed estimates and their bounds, 
#	and the graph name.
#	This function is used in the AnalyzeComplete function.
###########################################################################################
plotCIsSmooth = function(M, censortime, plotvalues, graphname, 
mainTitle = 'Hazard rate, smoothed', xlabel = 'Time', ylabel = 'Hazard estimate'){
	
	timePts = seq(censortime/2^M, censortime, by = censortime/2^M) - .5*censortime/2^M
	pdf(file = paste(graphname, '.pdf', sep = ''), width = 8, height = 6)
	plot(timePts, plotvalues[,2], main = mainTitle, xlab = xlabel, ylab = ylabel, type = 'l',
		 lty = 2, ylim = range(plotvalues))
	points(timePts, plotvalues[,1], type = 'l', lwd = 3)
	points(timePts, plotvalues[,3], type = 'l', lty = 2)
	dev.off()
}

###########################################################################################
# pbfxn is used to plot the results of the estimates across time.  It duplicates x and
#	y values appropriately so that it is easier to plot.
###########################################################################################
pbfxn = function(num.indic, unit.length, ydata){
	
	indx = rep(1:(num.indic+1), 2)
	indx = indx[order(indx)]
	
	xvals = (0:num.indic)[indx[-c(1,length(indx))]]*unit.length
	yvals = ydata[indx[-c(length(indx)-1, length(indx))]]
	
	outmat = as.data.frame(cbind(xvals, yvals))
	names(outmat) = c('x','y')
	
	return(outmat)
}

###########################################################################################
# findlogRatio is used to find the initial beta estimates for the non proportional 
#	hazards function.
###########################################################################################
#findlogRatio = function(Tivalues, trtGrpvalues, censvalues, Mval, censortimevalue){
	
#	fit0 = survfit(Surv(Tivalues[trtGrpvalues == 0], censvalues[trtGrpvalues == 0]) ~ 1)
#	fit1 = survfit(Surv(Tivalues[trtGrpvalues == 1], censvalues[trtGrpvalues == 1]) ~ 1)
	
#	hazRate0 = fit0$n.event/fit0$n.risk
#	hazRate1 = fit1$n.event/fit1$n.risk
	
#	timePts = (1:2^Mval)*censortimevalue/2^Mval
#	hazRate0 = sum(fit0$n.event[fit0$time <= timePts[1]]/fit0$n.risk[fit0$time <= timePts[1]])
#	hazRate1 = sum(fit1$n.event[fit1$time <= timePts[1]]/fit1$n.risk[fit1$time <= timePts[1]])
#	for(i in 2:length(timePts)){
#		hazRate0 = c(hazRate0, 
#					 sum(fit0$n.event[timePts[i-1]<fit0$time & fit0$time<=timePts[i]]/
#						 fit0$n.risk[timePts[i-1]<fit0$time & fit0$time<=timePts[i]]))
#		hazRate1 = c(hazRate1, 
#					 sum(fit1$n.event[timePts[i-1]<fit1$time & fit1$time<=timePts[i]]/
#						 fit1$n.risk[timePts[i-1]<fit1$time & fit1$time<=timePts[i]]))
#	}
	
#	lograts = log(hazRate1/hazRate0)
#	return(lograts)
#}

###########################################################################################
# prepX.categ is used to prep the categorical variables of the covariate matrix X,
#	converting factors to 0/1 variables.
###########################################################################################
prepX.categ = function(Xmatrix){
		
	# Loop through the covariate matrix, and identify all the variables that are factors, or there
	# are only two categories and they are not entered as 0/1.  factorIndex holds the columns that
	# contain factors, so once they have all been identified and converted to 0/1 columns, the
	# original columns can be deleted as they are not necessary for estimation.
	newX = NULL
	ncolXmat = 1
	nrowXmat = length(Xmatrix)
	if(!is.null(ncol(Xmatrix))){ 
		ncolXmat = ncol(Xmatrix)	
		nrowXmat = nrow(Xmatrix)
	}
	for(colctr in 1:ncolXmat){
		if(is.factor(as.data.frame(Xmatrix)[,colctr])){
			factorNames = unique(as.data.frame(Xmatrix)[,colctr])
			factorNames = factorNames[order(factorNames)]
			factorLength = length(factorNames)
			tempx = matrix(0, ncol = factorLength-1, nrow = nrowXmat)
			for(factorctr in 2:factorLength){
				tempx[which(as.data.frame(Xmatrix)[,colctr] == factorNames[factorctr]), factorctr-1] = 1
			}
			tempx = as.data.frame(tempx)
			tempxname = names(Xmatrix)[colctr]
			names(tempx) = paste(substr(tempxname, 1, nchar(tempxname)), factorNames[2:factorLength], sep = '.')
			if(colctr == 1){ newX = tempx
			} else { newX = cbind(newX, tempx)	}
		} else {
			newX = as.data.frame(cbind(newX, as.data.frame(Xmatrix)[,colctr]))
			names(newX) = c(names(newX)[-ncol(newX)], names(Xmatrix)[colctr])
		}
	}		
	return(as.data.frame(newX))
}

###########################################################################################
# prepX.contin is used to prep the continuous variables by standardizing them and
#	reducing their size if they are too large.  First, the ranges of each column in X
#	are calculated, and then they are multiplied by the lower and upper bounds.  If
#	any of the absolute values of these X*beta terms are larger than 75, the values
#	for that column are standardized (within column), and the new Xstd*beta ranges
#	are checked again.  If these absolute values are still larger than 75, then 
#	all the X values are divided by 5.  This continues until all |Xstd*beta| < 75.
###########################################################################################
prepX.contin = function(Xmatrix, blb, bub){
	
	# Find the parameters that are continuous by finding those that have 
	# more than two categories (factor variables are altered above, so 
	# this type of variable does not need to be taken in to account).	
	lengths = sapply(1:ncol(Xmatrix), function(x) length(table(Xmatrix[,x])[]))
	stdIndex = which(lengths > 2)
	Xmatrixmeans = Xmatrixsds = keepIndex = NULL
	innerctr = 0
	
	if(length(stdIndex) > 0){
		
		# Check the xbeta bounds, and standardize the Xs that do not fit in those bounds #
		ranges = sapply(stdIndex, function(x) range(Xmatrix[,x]))
		for(boundctr in 1:ncol(ranges)){
			
			absranges = max(abs(ranges[,boundctr]))*max(abs(c(blb[stdIndex[boundctr]], bub[stdIndex[boundctr]])))
			# If any of the Xbeta values are outside the range, then standardize them #
			if(absranges > 50){
		
				# Keep the index number
				keepIndex = c(keepIndex, stdIndex[boundctr])
				innerctr = innerctr + 1
			  
				# Standardize the X values
				Xmatrixmeans = c(Xmatrixmeans, mean(Xmatrix[,stdIndex[boundctr]]))
				Xmatrixsds = c(Xmatrixsds, sd(Xmatrix[,stdIndex[boundctr]]))
				Xmatrix[,stdIndex[boundctr]] = (Xmatrix[,stdIndex[boundctr]]-Xmatrixmeans[innerctr])/Xmatrixsds[innerctr]
		
				# Check the standardized X values to see if any of them are still too large
				absranges2 = max(abs(range(Xmatrix[,stdIndex[boundctr]])))*
					max(abs(c(blb[stdIndex[boundctr]], bub[stdIndex[boundctr]])))
				while(absranges2 > 50){
					Xmatrix[,stdIndex[boundctr]] = Xmatrix[,stdIndex[boundctr]]/5
					Xmatrixsds[innerctr] = Xmatrixsds[innerctr]*5
					absranges2 = max(abs(range(Xmatrix[,stdIndex[boundctr]])))*
						max(abs(c(blb[stdIndex[boundctr]], bub[stdIndex[boundctr]])))
				}
		    
			}
		}
	}
	
	return(list(X = Xmatrix, stdIndex = keepIndex, Xmeans = Xmatrixmeans, Xsds = Xmatrixsds))
}

###########################################################################################
# find.Rmp finds the Rmp values, given a vector of ds
###########################################################################################
find.Rmp = function(dvector){
	
	M = log(length(dvector))/log(2)
	loopM = M-1
	Rmp = NULL
	while(loopM >= 0){
		indexvec = (0:(2^(M-1-loopM)-1))*2^(loopM+1)+1
		Rmp = c(Rmp, 
				sapply(indexvec, function(x) 
					   sum(dvector[x:(x+2^loopM-1)])/sum(dvector[x:(x+2^(loopM+1)-1)])))
		loopM = loopM-1
	}
	return(Rmp)
}
		
###########################################################################################
# nph is used to get make a non-proportional hazards class for the formula in the library.
###########################################################################################
nph = function(x) x

###########################################################################################
# FindRoundVals is a function used to find the rounding values for a data set.
###########################################################################################
FindRoundVals = function(dataset){
	
	sdvalues = as.character(sapply(2:ncol(dataset), function(x) sd(dataset[,x])))
	round.index = (1:length(sdvalues))[which(sdvalues != 0)]
	
	# First split all characters to find the period.
	splitAt. = strsplit(sdvalues, '\\.')
	tempkeep = rep(0, length(sdvalues))
	for(ctr in round.index){	tempkeep[ctr] = splitAt.[[ctr]][2]	}
	
	# Find the actual rounding value.  It is either after the 'e-' in the number
	# or the number of zeroes after the decimal.
	round.values = rep(6, length(sdvalues))
	splitAte = strsplit(tempkeep, 'e-')
	for(ctr in round.index){
		if(length(splitAte[[ctr]]) > 1){
			temp.roundval = as.numeric(splitAte[[ctr]][2])
			if(temp.roundval > 5){	round.values[ctr] = temp.roundval + 1	}
		} else {
			temp.roundval = 0
			firstval = substr(tempkeep[ctr], 1, 1)
			while(firstval == '0'){
				temp.roundval = temp.roundval + 1
				tempkeep[ctr] = substr(tempkeep[ctr], 2, nchar(tempkeep[ctr]))
				firstval = substr(tempkeep[ctr], 1, 1)
			}
			if(temp.roundval > 6){	round.values[ctr] = temp.roundval + 2	}
		}
	}
	round.values[round.values > 16] = 16
	
	return(round.values)
}

###########################################################################################
# loglikePH calculates the -2*(log likelihood) value at each iteration of the MCMC chain
# for the PH model.
###########################################################################################
loglikePH = function(finaldata, numEstParams, censortime, Mval, numParams, delta, failBin, inBin, outfilename, X){

	nrow.data = nrow(finaldata)
	binwidth = censortime/2^Mval
    mat01 = calc_mat01(Mval)
    n = nrow(failBin)
	
	Rmpvec = as.data.frame(matrix(NA, ncol = sum(2^(1:Mval)), nrow = nrow.data))
	Rmpvec[,1:(sum(2^(1:Mval))/2)*2-1] = finaldata[,1:(2^Mval-1)+2^Mval+1+numParams+1]
	Rmpvec[,1:(sum(2^(1:Mval))/2)*2] = 1-finaldata[,1:(2^Mval-1)+2^Mval+1+numParams+1]
	
	logds = matrix(rep(log(finaldata$H00), each = 2^Mval), ncol = 2^Mval, nrow = nrow.data, byrow = TRUE) + 
					t(mat01%*%t(log(Rmpvec))) - log(binwidth)
	if(numParams > 0){	Xbeta = rowSums(as.matrix(X[rep(1:n, nrow.data),])*(as.matrix(finaldata[,1:numParams + 2^Mval+1])[rep(1:nrow.data, each = n),]))
	} else {    Xbeta = rep(0, n*nrow.data) }
	totalone = rep(delta, nrow.data)*(rowSums(failBin[rep(1:n, nrow.data),]*(logds[rep(1:nrow.data, each = n),]))+Xbeta) -
	rowSums(inBin[rep(1:n, nrow.data),]*(finaldata[,1:2^Mval+1][rep(1:nrow.data, each = n),]))*exp(Xbeta)
	
	logliks = -2*rowSums(matrix(totalone, ncol = n, byrow = TRUE))
	finaldata = as.data.frame(cbind(finaldata, logliks))
	names(finaldata) = c(names(finaldata)[-ncol(finaldata)], 'neg2loglike')
	write.table(finaldata, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
	
	AIC = 2*(numEstParams) - 2*min(logliks)
	BIC = -2*min(logliks) + (numEstParams)*log(n)
	DIC = .5*var(-2*logliks) + mean(-2*logliks)

	return(list(AIC = AIC, BIC = BIC, DIC = DIC))
}
	
###########################################################################################
# loglikeNPH calculates the -2*(log likelihood) value at each iteration of the MCMC chain
# for the NPH model.
###########################################################################################
loglikeNPH = function(finaldata, numEstParams, censortime, Mval, numParams, delta, failBin, inBin, outfilename, numHazards, Xfixed, numPHParams, indices){
	
	nrow.data = nrow(finaldata)
	binwidth = censortime/2^Mval
    n = nrow(failBin)
    mat01 = calc_mat01(Mval)
	Rmpvec = as.data.frame(matrix(NA, ncol = (2^Mval-1)*2*numHazards, nrow = nrow.data))
	Rmpvec[,1:((2^Mval-1)*2*numHazards/2)*2-1] = 
	finaldata[,1:((2^Mval-1)*numHazards) + 1+2^Mval*numHazards+numPHParams+2^Mval*(numHazards-1)+numHazards]
	Rmpvec[,1:((2^Mval-1)*2*numHazards/2)*2] = 
	1-finaldata[,1:((2^Mval-1)*numHazards) + 1+2^Mval*numHazards+numPHParams+2^Mval*(numHazards-1)+numHazards]
	logds = NULL
	for(hazctr in 1:numHazards){
		logds = cbind(logds, 
					  matrix(matrix(rep(log(finaldata[,1+numPHParams+2^Mval*(numHazards*2-1)+hazctr]), 
										each = 2^Mval), ncol = 2^Mval, nrow = nrow.data, byrow = TRUE)+ 
							 t(mat01%*%t(log(Rmpvec[,1:((2^Mval-1)*2)+(2^Mval-1)*2*(hazctr-1)]))) - 
							 log(binwidth), nrow = nrow.data))
	}
	if(numPHParams > 0){	
		Xbeta = rowSums(as.matrix(Xfixed[rep(1:n, nrow.data),])*(as.matrix(finaldata[,1:numPHParams+2^Mval*numHazards+1])[rep(1:nrow.data, each = n),]))
	} else { Xbeta = rep(0, n*nrow.data)	}
	hazBin = failBin[,rep(1:2^Mval, numHazards)]
	cumulBin = inBin[,rep(1:2^Mval, numHazards)]
	for(hazctr in 1:numHazards){
		hazBin[(1:n)[-indices[[hazctr]]], 1:2^Mval + (hazctr-1)*2^Mval] = 0
		cumulBin[(1:n)[-indices[[hazctr]]], 1:2^Mval + (hazctr-1)*2^Mval] = 0
	}
	totalone = rep(delta, nrow.data)*(rowSums(hazBin[rep(1:n, nrow.data),]*(logds[rep(1:nrow.data, each = n),]))+Xbeta) -
	rowSums(cumulBin[rep(1:n, nrow.data),]*(finaldata[,1:(2^Mval*numHazards)+1][rep(1:nrow.data, each = n),]))*exp(Xbeta)
	logliks = -2*rowSums(matrix(totalone, ncol = n, byrow = TRUE))
	finaldata = as.data.frame(cbind(finaldata, logliks))
	names(finaldata) = c(names(finaldata)[-ncol(finaldata)], 'neg2loglike')
	write.table(finaldata, paste(outfilename, '/MCMCchains.txt', sep = ''), row.names = FALSE)
	
	AIC = 2*(numEstParams) - 2*min(logliks)
	BIC = -2*min(logliks) + (numEstParams)*log(n)
	DIC = .5*var(-2*logliks) + mean(-2*logliks)
	
	return(list(AIC = AIC, BIC = BIC, DIC = DIC))
}
	

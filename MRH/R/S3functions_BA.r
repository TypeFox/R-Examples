###########################################################################################
# Contains the functions needed to: convert to an MRH object, summary, and plot for BA model.
###########################################################################################

callBAplot = function(x, censortime, main, xlab, ylab, plot.type, smooth.graph, smooth.df, 
combine.graphs, alpha.level, conf.int, log.ratio){

	firstgroup = substr(colnames(x)[2], 4, 4)
	dnames = colnames(x)[which(substr(colnames(x), 1, 1) == 'd')]
	numbins = length(which(unlist(strsplit(dnames, '_'))[1:length(dnames)*2] == firstgroup))
	M = log(numbins)/log(2)
	hazgroups = unique(unlist(strsplit(dnames, '_'))[1:length(dnames)*2])
	numhazards = length(hazgroups)
	
	# Get the parameter estimates #
	if(plot.type == 'r'){
		if(log.ratio != TRUE){
			start.betas = substr(colnames(x), 1, 4)
			end.betas = substr(colnames(x), nchar(colnames(x))-3, nchar(colnames(x)))
			beta.nphindex = 1:(2^M*(numhazards-1)) + which(start.betas == 'beta' & end.betas == 'bin1')[1] - 1
			x[,beta.nphindex] = exp(x[,beta.nphindex])
		}
        dests.temp = summary.MRH(x, alpha.level = alpha.level, maxStudyTime = censortime)$beta
        end.betas = substr(row.names(dests.temp), nchar(row.names(dests.temp))-3, nchar(row.names(dests.temp)))
        beta.nph.index = 1:(2^M*(numhazards-1))+which(end.betas == 'bin1')[1]-1
        dests.temp = dests.temp[beta.nph.index,]
	} else if(!missing(censortime)){
		dests.temp = t(CalcFunction(x, function.type = plot.type, alpha.level = alpha.level, maxStudyTime = censortime)[[1]])[,c(2,1,3)]
	} else {	dests.temp = t(CalcFunction(x, function.type = plot.type, alpha.level = alpha.level)[[1]])[,c(2,1,3)]	}	
	if(plot.type %in% c('H','S')){	dests.temp = dests.temp[-((1:numhazards-1)*(2^M+1)+1),]	}
	if(smooth.graph == TRUE){	
		if(is.null(smooth.df)){	smooth.df = 2^M/2	}
		smooth.timepts = 1:2^M*censortime/2^M-.5*censortime/2^M
		dests = NULL
		for(ctr in 1:(numhazards-1)){
			dests = rbind(dests, 
							cbind(smooth.spline(dests.temp[1:2^M+2^M*(ctr-1),1]~smooth.timepts, df = smooth.df)$y,
								  smooth.spline(dests.temp[1:2^M+2^M*(ctr-1),2]~smooth.timepts, df = smooth.df)$y,
								  smooth.spline(dests.temp[1:2^M+2^M*(ctr-1),3]~smooth.timepts, df = smooth.df)$y))
		}
		if(plot.type != 'r'){
			ctr = numhazards
			dests = rbind(dests, cbind(smooth.spline(dests.temp[1:2^M+2^M*(ctr-1),1]~smooth.timepts, df = smooth.df)$y,
									   smooth.spline(dests.temp[1:2^M+2^M*(ctr-1),2]~smooth.timepts, df = smooth.df)$y,
									   smooth.spline(dests.temp[1:2^M+2^M*(ctr-1),3]~smooth.timepts, df = smooth.df)$y))
		} 		
		plot.timepts = smooth.timepts
		max.index = 2^M
	} else {	
		plot.timepts = pbfxn(2^M, censortime/2^M, rep(NA, 2^M))$x	
		dests = NULL
		for(ctr in 1:(numhazards-1)){
			dests = rbind(dests, 
						  cbind(pbfxn(2^M, censortime/2^M, dests.temp[1:2^M+2^M*(ctr-1),1])$y, 
								pbfxn(2^M, censortime/2^M, dests.temp[1:2^M+2^M*(ctr-1),2])$y,
								pbfxn(2^M, censortime/2^M, dests.temp[1:2^M+2^M*(ctr-1),3])$y))
		}
		if(plot.type != 'r'){
			ctr = numhazards
			dests = rbind(dests, 
						  cbind(pbfxn(2^M, censortime/2^M, dests.temp[1:2^M+2^M*(ctr-1),1])$y, 
								pbfxn(2^M, censortime/2^M, dests.temp[1:2^M+2^M*(ctr-1),2])$y,
								pbfxn(2^M, censortime/2^M, dests.temp[1:2^M+2^M*(ctr-1),3])$y))
		}
		max.index = 2^M*2
	}
	
	# Make the y labels and get the appropriate x and y values.
	ylabel = 'Hazard rate'
	if('H' == plot.type){	
		ylabel = 'Cumulative hazard'
	} else if('S' == plot.type){	
		ylabel = 'Survival function'
	} else if('r' == plot.type){
		if(log.ratio == TRUE){
			ylabel = 'Log hazard ratio'	
		} else {
			ylabel = 'Hazard ratio'
		}
	}

	timeptslabel = plot.timepts
	
	###########################################################################################################
	#								SEPARATE ESTIMATES
	###########################################################################################################
	ylimit = range(dests)
	if(conf.int != TRUE){ ylimit = range(dests[,1])	}
	if(plot.type == 'S'){ ylimit[1] = 0	}
	if(combine.graphs != TRUE){
		if('r' != plot.type){
			plotnums = c(ceiling(numhazards/2), 2)
			numplots = numhazards
			mainnames = paste(hazgroups, 'group')
			mtextname = 'Hazard Rates by Group'
		} else { 
			if(numhazards > 2){	plotnums = c(ceiling((numhazards-1)/2), 2)	
			} else {	plotnums = c(1,1)	}
			numplots = numhazards-1
			mainnames = paste('Group', hazgroups[-1])
			if(log.ratio == TRUE){
				mtextname = paste('Log hazard Ratios: Comparison to Group', hazgroups[1])
			} else {
				mtextname = paste('Hazard Ratios: Comparison to Group', hazgroups[1])
			}
		}
		par(mfrow = plotnums, mai = c(.35, .1, .25, 0), oma = c(0, 2, 2, 1), tck = -.02)
		for(hazCtr in 1:numplots){
			plot(plot.timepts, dests[1:max.index+max.index*(hazCtr-1),1], type = 'l', lwd = 3, ylab = '', xlab = '', 
				 ylim = ylimit, main = mainnames[hazCtr], axes = FALSE)
			if(plot.type == 'r'){	
				if(log.ratio == TRUE){	abline(h = 0, col = 'grey')	} else {	abline(h = 1, col = 'grey')	}
			}
			if(conf.int == TRUE){
				points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1),2], type = 'l', lty = 2)
				points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1),3], type = 'l', lty = 2)
			}
			box()
			axis(1, at = timeptslabel, labels = round(timeptslabel, 1), cex.axis = .75, padj = -2)
			mtext('Time', side = 1, padj = 1.5, cex = .8)
			if(hazCtr %% 2 != 0){	axis(2, padj = 1.25)	}
		}
		mtext(mtextname, outer = TRUE, cex = 1.4)
	###########################################################################################################
	#								COMBINED ESTIMATES
	###########################################################################################################
	} else { 
		plot(plot.timepts, dests[1:max.index,1], type = 'l', lwd = 2,
			 xlab = 'Time', ylab = ylabel, main = ylabel, ylim = ylimit)
		if(plot.type == 'r'){
			if(log.ratio == TRUE){	abline(h = 0, col = 'grey')	} else {	abline(h = 1, col = 'grey')	}
		}
		if(conf.int == TRUE){
			points(plot.timepts, dests[1:max.index,2], type = 'l', lty = 2)
			points(plot.timepts, dests[1:max.index,3], type = 'l', lty = 2)
		}
		if((plot.type != 'r') | (plot.type == 'r' & numhazards > 2)){
		   for(hazCtr in 2:(numhazards-1)){
				points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1), 1], type = 'l', lwd = 2, col = hazCtr)
				if(conf.int == TRUE){
					points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1), 2], type = 'l', lty = 2, col = hazCtr)
					points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1), 3], type = 'l', lty = 2, col = hazCtr)
				}
		   }
		}
		if(plot.type != 'r'){
			hazCtr = numhazards
			points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1), 1], type = 'l', lwd = 2, col = hazCtr)
			if(conf.int == TRUE){
				points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1), 2], type = 'l', lty = 2, col = hazCtr)
				points(plot.timepts, dests[1:max.index+max.index*(hazCtr-1), 3], type = 'l', lty = 2, col = hazCtr)
			}
			legend(x = "topright", paste('group', hazgroups), fill = 1:numhazards, cex = 1.2)
		} else { 
			legend(x = "topright", paste('group', hazgroups[-1]), fill = 2:numhazards-1, cex = 1.2)
		}
	}
}

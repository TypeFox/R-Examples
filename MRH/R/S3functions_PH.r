###########################################################################################
# Contains the functions needed to: convert to an MRH object, summary, and plot for PH model.
###########################################################################################

callPHplot = function(x, censortime, main, xlab, ylab, plot.type, 
smooth.graph, smooth.df, combine.graphs, conf.int, alpha.level){
	
	# Calculate M 
	M = log(length(which(substr(colnames(x), 1, 1) == 'd')))/log(2)
	# Cannot smooth with 2 bins, so set smooth = F if user has selected this option 
	# (the graph will still appear smooth to the user)
	if(M == 1 & smooth.graph == TRUE){	smooth.graph = FALSE	}

	# Get the parameter estimates #
	results = t(CalcFunction(x, function.type = plot.type, alpha.level = alpha.level, maxStudyTime = censortime)[[1]])
	
	if(plot.type == 'H'){	ylab = 'Cumulative hazard'
	} else if(plot.type == 'S'){	ylab = 'Survival Function'	}
	
	###### Plot #####
	if(smooth.graph == FALSE){
		if(plot.type == 'h'){	xvals = pbfxn(2^M, censortime/2^M, as.numeric(results[,2]))
		} else {
			xvals = pbfxn(2^M, censortime/2^M, results[-1,2])
			if(plot.type == 'H'){	xvals = rbind(c(0, 0), xvals)
			} else { xvals = rbind(c(0, 1), xvals)	}
		}
		ylimit = range(results)
		if(conf.int == FALSE){	ylimit = range(results[,2])	}
		plot(xvals$x, xvals$y, main = main, xlab = xlab, ylab = ylab, type = 'l', 
			 ylim = ylimit, lwd = 3)
		if(conf.int == TRUE){
			if(plot.type == 'h'){
				points(xvals$x, pbfxn(2^M, censortime/2^M, results[,1])$y, type = 'l', lty = 2)
				points(xvals$x, pbfxn(2^M, censortime/2^M, results[,3])$y, type = 'l', lty = 2)
			} else if(plot.type == 'H'){
				points(xvals$x, c(0, pbfxn(2^M, censortime/2^M, results[-1,1])$y), type = 'l', lty = 2)
				points(xvals$x, c(0, pbfxn(2^M, censortime/2^M, results[-1,3])$y), type = 'l', lty = 2)
			} else {
				points(xvals$x, c(1, pbfxn(2^M, censortime/2^M, results[-1,1])$y), type = 'l', lty = 2)
				points(xvals$x, c(1, pbfxn(2^M, censortime/2^M, results[-1,3])$y), type = 'l', lty = 2)
			}				
		}
	} else {	
		# Calculate and plot the smooth ds over time 
		if(is.null(smooth.df)){	smooth.df = 2^M/2	}
		timePts = seq(censortime/2^M, censortime, by = censortime/2^M)-.5*censortime/2^M
		if(plot.type != 'h'){	timePts = c(0, timePts)	}
		smoothdInts = smooth.spline(results[,1]~timePts, df = smooth.df)$y
		smoothdInts = cbind(smoothdInts, smooth.spline(results[,2]~timePts, df = smooth.df)$y)
		smoothdInts = cbind(smoothdInts, smooth.spline(results[,3]~timePts, df = smooth.df)$y)
		ylimit = range(smoothdInts)
		if(conf.int == FALSE){	ylimit = range(smoothdInts[,2])	}
		plot(timePts, smoothdInts[,2], main = main, xlab = xlab, ylab = ylab, type = 'l', 
			 ylim = ylimit, lwd = 3)
		if(conf.int == TRUE){
			points(timePts, smoothdInts[,1], type = 'l', lty = 2)
			points(timePts, smoothdInts[,3], type = 'l', lty = 2)
		}
	}
}	

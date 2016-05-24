#
# tool functions
#

# returns a set of colors used by ggplot2, they are pretty good
.ggplotcolors <- function(n=6, h=c(0, 360) +15){
	if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
	hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

# draw the variance shading around mean
.drawShading <- function(x, mean, std, color, noisestd=NULL, 
								 minalpha=0.05, maxalpha=0.55, sigma=2, plotgradient=TRUE) {
	r <- col2rgb(color)[1]/255
	g <- col2rgb(color)[2]/255
	bl <- col2rgb(color)[3]/255
	
	pupper95 = mean + sigma*std
	plower95 = mean - sigma*std
	
	# draw posterior 95% CI
	lines(x=c(x,rev(x)), y=c(plower95,rev(pupper95)), col=rgb(r,g,bl,maxalpha))
	# also draw noise curves
	if (!is.null(noisestd)) {
		npupper95 = mean + sigma*(sqrt(std^2 + noisestd^2))
		nplower95 = mean - sigma*(sqrt(std^2 + noisestd^2))
		lines(x=c(x,rev(x)), y=c(nplower95,rev(npupper95)), col=rgb(r,g,bl,1), lty=2, lwd=1.3)
	}
	
	# draw shading color by variance
	if (plotgradient) {
		for (i in 1:length(x)) {
			if (i == 1) {
				xprev <- x[1]
				ylprev <- plower95[1]
				yuprev <- pupper95[1]
				xnext <- mean(c(x[i],x[i+1]))
				ylnext <- mean(c(plower95[i], plower95[i+1]))
				yunext <- mean(c(pupper95[i], pupper95[i+1]))
			} else if (i == length(x)) {
				xprev <- mean(c(x[i],x[i-1]))
				ylprev <- mean(c(plower95[i], plower95[i-1]))
				yuprev <- mean(c(pupper95[i], pupper95[i-1]))
				xnext <- x[i]
				yunext <- pupper95[i]
				ylnext <- plower95[i]
			} else {
				# mean time point between current and previous
				xprev <- mean(c(x[i],x[i-1]))
				xnext <- mean(c(x[i],x[i+1]))
				ylprev <- mean(c(plower95[i], plower95[i-1]))
				ylnext <- mean(c(plower95[i], plower95[i+1]))
				yuprev <- mean(c(pupper95[i], pupper95[i-1]))
				yunext <- mean(c(pupper95[i], pupper95[i+1]))
			}    
			
			polygon(x= c(xprev, xnext, xnext, xprev), 
					  y= c(ylprev, ylnext, yunext, yuprev), 
					  density=NA, col= rgb(r,g,bl, max(minalpha,maxalpha-1.1*std[i])), border=NA)
		}
	} else {  # draw simple polygon
		polygon(x= c(x,rev(x)), 
				  y= c(mean-2*std,rev(mean)+2*rev(std)), 
				  density=NA, col= rgb(r,g,bl, (minalpha+maxalpha)/2), border=NA)
	}
}

#' @import MASS
.drawSamples <- function(timepoints, mean, cov, N=3) {
	sample = mvrnorm(N, mean, cov)
	for (i in 1:N) {
		lines(x= timepoints, y= sample[i,], col = "grey", lw=1)
		points(timepoints, sample[i,], col='black', pch=20, cex=0.6)
	}
}

# add lengthscale on right y-axis
.plotls = function(obj) {
	par(new=T)
	plot(obj$targets$x, .bexp(obj$targets$x, obj$params$l, obj$params$lmin, obj$params$c), axes=F, ylab="", ylim=c(0,6),
		  t='l', lty=1, lwd=1.2, xaxs='i', col='blue', xlab='')
	axis(4, at=c(1,3,5), col.axis='blue', mgp=c(3,0.5,0))
	mtext('ell(t)', side=4, line=1.5, col='blue', cex=1.2)
}


#' @keywords internal
#' @title Plot a simple gaussian process
#' @description Plots a gaussian process. Several boolean parameters for modifying the 
#'  plot. By default plots the data, posterior mean and 95\% interval. 
#' @param x the gp-object
#' @param y placeholder variable
#' @param plotdata plot the data (default)
#' @param plotmean plot the GP mean (default)
#' @param plotcov plot the GP covariances (default)
#' @param plotnoise plot the observational noise (default)
#' @param samples plot N samples from the GP
#' @param sigma variance level to plot
#' @param title plot title
#' @param legend plot legend
#' @param plotgradient use gradient graphics
#' @param ... ...
plot.gpsimple = function(x, y=NULL, plotdata=TRUE, plotmean=TRUE, plotcov=TRUE, plotnoise=FALSE, 
								 samples=0, sigma=2, title=NULL, legend=FALSE, plotgradient=TRUE,...) {
	obj=x
	xmin <- min(obj$x)
	xmax <- max(obj$x)
	
	if (!hasArg("ylim") & plotnoise == TRUE) {
		ymax = max(sigmaupper(obj,T,sigma))
		ymin = min(sigmalower(obj,T,sigma))
		ylim = c(ymin-0.2,ymax+0.2)
	} else if (!hasArg("ylim") & plotnoise == FALSE) {
		ymax = max(sigmaupper(obj,F,sigma))
		ymin = min(sigmalower(obj,F,sigma))
		ylim = c(ymin-0.2,ymax+0.2)
	} else {
		ylim = list(...)$ylim
	}
	
	if (!hasArg("xlab")) {  xlab = 'Time'  } else { xlab = list(...)$xlab }
	if (!hasArg("ylab")) {	ylab = 'Value' } else { ylab = list(...)$ylab }
	
	colors <- .ggplotcolors(n=3)[1]
	plot(NA, xlim=c(xmin,xmax), ylim=ylim, xlab='', ylab='', xaxs='i', main=title, mgp=c(3,0.7,0))
	mtext(xlab, side=1, line=2, cex=1.2)
	mtext(ylab, side=2, line=1.8, cex=1.2)	
	
	if (plotcov && plotnoise) {
		.drawShading(obj$x, obj$mean, sqrt(diag(obj$cov)), colors, obj$noisestd, plotgradient=plotgradient, sigma=sigma)
	} else if (plotcov) {
		.drawShading(obj$x, obj$mean, sqrt(diag(obj$cov)), colors, plotgradient=plotgradient, sigma=sigma)
	}
	if (plotmean)       { lines(obj$x, obj$mean, col=colors, lwd=3) }
	if (plotdata)       { points(obj$x.obs, obj$y.obs, bg=colors, pch=21) }
	if (samples > 0)    { .drawSamples(obj$x, obj$mean, obj$cov, samples) }
	if (legend) { 
		legend('topright', c('Mean', paste('sigma=',sigma,' CI',sep=''),paste('sigma=',sigma,' noisy CI',sep='')), lty=c(1,1,2), lwd=c(3,1,2), col=c('red','red','grey'))
	}
}

#' @export
#' @title Plot a gaussian process
#' @description Plots a gaussian process. Several boolean parameters for modifying the 
#'  plot. By default plots the data, posterior mean and 95\% interval. 
#' @param x the gp-object
#' @param y placeholder variable
#' @param plotdata plot the data (default)
#' @param plotmean plot the GP mean (default)
#' @param plotcov plot the GP covariances (default)
#' @param plotnoise plot the observational noise (default)
#' @param samples plot N samples from the GP
#' @param sigma variance level to plot
#' @param title plot title
#' @param legend plot legend
#' @param plotgradient use gradient graphics
#' @param plotls plot lengthscale function
#' @param ... ...
#' @examples
#' # read toy data
#' data(toydata)
#'
#' \dontrun{can take several minutes
#' # perform one-sample regression
#' res = gpr2sample(toydata$ctrl$x, toydata$ctrl$y, seq(0,22,0.1))
#'	
#' # pre-computed model for toydata
#' data(toygps)
#' res = toygps$ctrlmodel
#'	
#' # basic plot with data, estimated mean and 95\%
#' plot(res)
#'	
#' # don't plot the data, plot some samples drawn from the learned gp
#' plot(res, plotdata=FALSE, samples=3)}
plot.gp <- function(x, y=NULL, plotdata=TRUE, plotmean=TRUE, plotcov=TRUE, plotnoise=FALSE, 
						  samples=0, sigma=2, title=NULL, legend=FALSE, plotgradient=TRUE, plotls=FALSE, ...) {
	obj=x
	xmin <- min(obj$x)
	xmax <- max(obj$x)
	
	if (!hasArg("ylim") & plotnoise == TRUE) {
		ymax = max(sigmaupper(obj,T,sigma))
		ymin = min(sigmalower(obj,T,sigma))
		ylim = c(ymin-0.2,ymax+0.2)
	} else if (!hasArg("ylim") & plotnoise == FALSE) {
		ymax = max(sigmaupper(obj,F,sigma))
		ymin = min(sigmalower(obj,F,sigma))
		ylim = c(ymin-0.2,ymax+0.2)
	} else {
		ylim = list(...)$ylim
	}
	
	if (!hasArg("xlab")) {  xlab = 'Time'  } else { xlab = list(...)$xlab }
	if (!hasArg("ylab")) {	ylab = 'Value' } else { ylab = list(...)$ylab }
	
	colors <- .ggplotcolors(n=3)[1]
	plot(NA, xlim=c(xmin,xmax), ylim=ylim, xlab='', ylab='', xaxs='i', main=title, mgp=c(3,0.7,0))
	mtext(xlab, side=1, line=2, cex=1.2)
	mtext(ylab, side=2, line=1.8, cex=1.2)	
	
	if (plotcov && plotnoise) {
		.drawShading(obj$targets$x, obj$targets$pmean, obj$targets$pstd, colors, obj$targets$noisestd, plotgradient=plotgradient, sigma=sigma)
	} else if (plotcov) {
		.drawShading(obj$targets$x, obj$targets$pmean, obj$targets$pstd, colors, plotgradient=plotgradient, sigma=sigma)
	}
	if (plotmean)       { lines(obj$targets$x, obj$targets$pmean, col=colors, lwd=3) }
	if (plotdata)       { points(obj$x, obj$y, bg=colors, pch=21) }
	if (samples > 0)    { .drawSamples(obj$targets$x, obj$targets$pmean, obj$cov, samples) }
	if (legend) { 
		legend('topright', c('Mean', paste('sigma=',sigma,' CI',sep=''),paste('sigma=',sigma,' noisy CI',sep='')), lty=c(1,1,2), lwd=c(3,1,2), col=c('red','red','grey'))
	}
	
	# draw also the lengthscale function
	if (plotls) { .plotls(obj) }
}

#' @export
#' @title Plots several GP's simultaneously
#' @description Plots the GP's corresponding to the control and case data, 
#'  as well as the null model. Visualizes the log likelihood ratios between 
#'  the null and individual models. Several boolean parameters for modifying
#'  the plot. By default plots the data, posterior mean and 95\% interval for 
#'  CASE and CONTROL.
#' @details The threshold \code{thr} is the logarithmic likelihood ratio between
#'  null and control+case models. The default value 1 hence corresponds to 
#'  a likelihood ratio of 2.72.
#' @param x the \code{gppack}-object
#' @param y placeholder variable
#' @param plotdata plot the data (default)
#' @param plotmeans plot the GP mean (default)
#' @param plotcovs plot the GP covariances (default)
#' @param plotnoises plot the observational noise (default)
#' @param plotnull plots also the null model
#' @param plotratios plots the ratios, choices are \code{emll}, \code{mll}, \code{npc}, \code{pc}
#' @param thr ratio threshold
#' @param samples plot N samples from the GP
#' @param sigma variance level to plot
#' @param title plot title
#' @param legend plot legend
#' @param plotgradient use gradient graphics
#' @param ... ...
#' @examples
#' # read toy data
#' data(toydata)
#'	
#' \dontrun{can take several minutes
#' # perform two-sample regression
#' res = gpr2sample(toydata$ctrl$x, toydata$ctrl$y, toydata$case$x, toydata$case$y, seq(0,22,0.1))
#'	
#' # pre-computed model for toydata
#' data(toygps)
#' res = toygps
#'	
#' # basic plot
#' plot(res)
#'	
#' # plot also the null model, don't plot data, means or noise
#' plot(res, plotnull=TRUE, plotdata=FALSE, plotmeans=FALSE, plotnoise=FALSE)}
plot.gppack <- function(x, y=NULL, plotdata=TRUE, plotmeans=TRUE, plotcovs=TRUE, plotnoises=FALSE,
								plotnull=FALSE, plotratios='emll', thr=1, 
								samples=0, sigma=2, title=NULL, legend=FALSE, plotgradient=TRUE, ...) {
	obj=x
	x.targets = obj$ctrlmodel$targets$x
	K <- length(x.targets) # timepoint count
	xmin <- min(x.targets)
	xmax <- max(x.targets)
	
	if (!hasArg("xlab")) {  xlab = 'Time'  } else { xlab = list(...)$xlab }
	if (!hasArg("ylab")) {	ylab = 'Value' } else { ylab = list(...)$ylab }
	
	colors <- .ggplotcolors(n=3)
	colors = colors[c(3,2,1)]
	
	if (plotcovs & plotnoises) {
		maxcase = max(sigmaupper(obj$casemodel,T,sigma))
		maxctrl = max(sigmaupper(obj$ctrlmodel,T,sigma))
		mincase = min(sigmalower(obj$casemodel,T,sigma))
		minctrl = min(sigmalower(obj$ctrlmodel,T,sigma))
		ymax = max(c(maxctrl,maxcase))
		ymin = min(c(minctrl,mincase))
	} else if (plotcovs) {
		ymax <- max(max(c(sigmaupper(obj$ctrlmodel,F,sigma), sigmaupper(obj$casemodel,F,sigma))))
		ymin <- min(min(c(sigmalower(obj$casemodel,F,sigma), sigmalower(obj$casemodel,F,sigma))))
	} else {
		ymax = max(max(c(obj$ctrlmodel$y, obj$casemodel$y)))
		ymin = min(min(c(obj$ctrlmodel$y, obj$casemodel$y)))
	}
	
	################### upper
	if (plotratios == "emll" || plotratios == 'npc' || plotratios == 'pc' || plotratios == 'mll') {
		# Use base-graphics
		layout(matrix(c(1,2), 2, 1, byrow=TRUE), heights=c(0.15,0.85), widths=1)
		
		par(mar=c(0,4,1,1), mgp=c(3,0.8,0))
		plot(-1, xlim=c(xmin, xmax), ylim=c(0, max(2,2*thr)), ylab="", xaxt="n", 
			  yaxt="n",xaxs="i", yaxs="i", main=title, cex.main=1.2)
		
		mtext(text="log\nratio", side=2, line=1.4, cex=1.4)
		axis(2, at=c(0, thr, max(2, 2*thr)), las=2, cex.axis=.9, tck=0.06)
		axis(4, tck=0.06, at=c(1, thr, 3), labels=NA)
		axis(1, tck=0.06, labels=NA)
		axis(3, tck=0.06, labels=NA)
		
		lines(x.targets, obj$ratios[[plotratios]], lwd=1.5)
		polygon(c(x.targets,rev(x.targets)), 
				  c(rep(0,K), rev(obj$ratios[[plotratios]])), 
				  density=NA, border=NA, col=colors[3])  
		
		if (!is.null(thr)) {
			lines(x=x.targets, y=rep(thr,K), col="grey", lwd=1)
		}
	}
	
	
	if (!hasArg("ylim")) {	ylim = c(ymin-0.2,ymax+0.2)  } else {	ylim = list(...)$ylim }
	
	################### lower
	par(mar=c(5,4,0,1))
	plot(NA, xlim=c(xmin,xmax), ylim=ylim, xlab='', ylab='', xaxt="n", yaxt="n", xaxs="i", yaxs="i")
	
	mtext(text=xlab, side=1, line=2.5, cex=1.4)
	mtext(text=ylab, side=2, line=2.5, cex=1.4)
	axis(2, las=2, cex.axis=.9, tck=0.015)
	axis(1, tck=0.015)
	axis(4, labels=NA, tck=0.015)
	axis(3, tck=0.015, labels=NA)
	
	if (legend==TRUE) {
		if (plotnoises && plotnull) {
			legend('topright', c('Control posterior', 'Case posterior', 'Shared posterior', '  noisy posterior'), pch=c(20,20,NA,NA), col=c(colors[1],colors[3],colors[2],'grey45'), lty=c(1,1,1,2))
		} else if (!plotnoises && plotnull) {
			legend('topright', c('Control posterior', 'Case posterior', 'Shared posterior'), pch=c(20,20,NA), col=c(colors[1],colors[3],colors[2]), lty=c(1,1,1))
		} else if (plotnoises && !plotnull) {
			legend('topright', c('Control posterior', 'Case posterior', '  noisy posterior'), pch=c(20,20,NA), col=c(colors[1],colors[3],'gray45'), lty=c(1,1,2))
		} else {
			legend('topright', c('Control posterior', 'Case posterior'), pch=c(20,20), col=c(colors[1],colors[3]), lty=c(1,1))
		}
	}
	
	if (!is.null(thr) && (plotratios == "emll" || plotratios == 'npc' || plotratios == 'pc' || plotratios == 'mll')) {
		i <- 1
		while(i < K) {
			if (obj$ratios[[plotratios]][i]>thr) {
				j <- 1
				while(obj$ratios[[plotratios]][i+j]>thr && i+j < K-1){
					j <- j+1
				}
				xprev <- mean(c(x.targets[i],x.targets[i-1]))
				xnext <- mean(c(x.targets[i+j],x.targets[i+j+1]))
				
				polygon(x=c(xprev, xnext, xnext, xprev), 
						  y=c(-100,-100,100,100), 
						  density=NA, border=NA, col=rgb(0.0,0.0,0.0,0.10))
				i <- i+j+1
			}
			else { i <- i+1 }
		}
	}
	
	if (plotcovs) {
		if (plotnoises) {
			.drawShading(x.targets, obj$ctrlmodel$targets$pmean, obj$ctrlmodel$targets$pstd, colors[1], obj$ctrlmodel$targets$noisestd, plotgradient=plotgradient)
			.drawShading(x.targets, obj$casemodel$targets$pmean, obj$casemodel$targets$pstd, colors[3], obj$casemodel$targets$noisestd, plotgradient=plotgradient)
			if (plotnull) {
				.drawShading(x.targets, obj$nullmodel$targets$pmean, obj$nullmodel$targets$pstd, colors[2], obj$nullmodel$targets$noisestd, plotgradient=plotgradient)
			}
		} else {
			.drawShading(x.targets, obj$ctrlmodel$targets$pmean, obj$ctrlmodel$targets$pstd, colors[1], plotgradient=plotgradient)
			.drawShading(x.targets, obj$casemodel$targets$pmean, obj$casemodel$targets$pstd, colors[3], plotgradient=plotgradient)
			if (plotnull) {
				.drawShading(x.targets, obj$nullmodel$targets$pmean, obj$nullmodel$targets$pstd, colors[2], plotgradient=plotgradient)
			}
		}
	}
	
	if (plotdata) {
		points(obj$ctrlmodel$x, obj$ctrlmodel$y, bg=colors[1], pch=21, cex=1.4)
		points(obj$casemodel$x, obj$casemodel$y, bg=colors[3], pch=21, cex=1.4)
	}
	
	if (plotmeans) {
		lines(x.targets, obj$ctrlmodel$targets$pmean, col=colors[1], lwd=3.5)
		lines(x.targets, obj$casemodel$targets$pmean, col=colors[3], lwd=3.5)
		if (plotnull) {
			lines(x.targets, obj$nullmodel$targets$pmean, col=colors[2], lwd=3.5)
		}
	}
}





# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # draw error bars
# .errorBars = function(x, y, noise, sigma=2, bar=0.1) {
#   segments(x,y-sigma*noise,x,y+sigma*noise)
#   segments(x-bar,y-sigma*noise,x+bar,y-sigma*noise)
#   segments(x-bar,y+sigma*noise,x+bar,y+sigma*noise)
# }
# 
# 
# # Draw a colored line(or area), Red/Green/Black(or transparent) based on (mean +/- 2*std)'s value
# .drawRGB <- function(timepoints, mean, std, ylimits, lwd=3, shift=0){
#   # Shift value is used to draw multiple curves in a single plot by shifting it's y value.
#   color <- ifelse( (mean+2*std)<(1+shift),'forestgreen',ifelse( (mean-2*std)>(1+shift),'red','black'))
#   # Colored lines
#   #   s <- seq(length(timepoints)-1)
#   # 	segments(timepoints[s], mean[s],timepoints[s+1], mean[s+1], col= color, lwd=lwd)
#   #   segments(timepoints[s], 1+shift,timepoints[s+1], 1+shift  , col= color, lwd=lwd/2)
#   # Colored areas
#   i <- 1
#   while(i < (length(timepoints)-1)){
#     actualcol <- color[i]
#     j <- 1
#     while(i+j <= (length(timepoints)-1) && color[i+j] == actualcol){
#       j <- j+1
#     }
#     if(actualcol!="black"){
#       polygon(x=c(timepoints[i], timepoints[i+j], timepoints[i+j], timepoints[i]), 
#               y=c(ylimits[1]+shift,ylimits[1]+shift,ylimits[2]+shift,ylimits[2]+shift), 
#               density=NA, border=NA, 
#               col=rgb(col2rgb(actualcol)[1]/255,col2rgb(actualcol)[2]/255,col2rgb(actualcol)[3]/255,0.25))
#     }
#     i <- i+j+1
#   }
# }
# 
# # Adds \n in a text
# .wrap.it <- function(x, len)
# { 
#   sapply(x, function(y) paste(strwrap(y, len), 
#                               collapse = "\n"), 
#          USE.NAMES = FALSE)
# }
# 
# ############################################################
# #################### PLOTTING FUNCTIONS ####################
# ############################################################
# 
# 
# 
# 
# # Create a set of gp plots. Either the whole data in a PDF (several pages) or just the first page in R.
# # The "names" argument require the full csv file since it shows the prot name, id and gene name.
# plot.gpset <- function(x, names=NULL, nrow=6, ncol=6, file=TRUE, filename="gpsetPlot"
# 					   , plotdata=FALSE, plotnoise=FALSE, colored=FALSE, uniformvar=FALSE, lwd=3, A4=TRUE){
# 	#  nrow=max(1,min(50,nrow))
# 	ncol=max(1,min(50,ncol))
# 	ext = ".pdf"
# 	filename = paste0(filename,ext)
# 	obj = x
# 	objSize = length(obj)
# 	if(!colored){
# 		color <- .ggplotcolors(n=3)[1]
# 		colorLine <- color
# 	} else{
# 		color <- "gray"
# 		colorLine <- "black"
# 	}
# 	txt <- 30 # White area where text is displayed
# 	ysize <- 50 # Height of the plot
# 	minSpace <- 0.4 # Minimum space between 2 rows of curves
# 	if(A4){
# 		lwd <- 1.5
# 		ysize <- 30
# 		# Settings without protNames
# 		ncol <- 5
# 		# Settings with protNames
# 		#     ncol <- 4
# 		#     txt <- 70
# 		#     minSpace <- 0.7
# 	}
# 	count <- 1
# 	# Order the data with the biggest plots first
# 	ord<-c()
# 	for(i in 1:length(obj)){
# 		ord<-c(ord,max((1-min(obj[[i]]$targets$plower95)),(max(obj[[i]]$targets$pupper95)-1)))
# 	}
# 	obj <- obj[order(ord,decreasing=T)]
# 	if(!is.null(names)){
#     if(nrow(names)==length(ord)){
#   		geneName <- names[order(ord,decreasing=T),"symbol"]
#   		protAcc <- names[order(ord,decreasing=T),"uniprot"]
#   		protName <- names[order(ord,decreasing=T),"prot_name"]
#   	}
# 	}
# 	ord <- sort(ord,decreasing=T)
# 	
# 	# Start of the plotting.
# 	if(file){
# 		if(A4){
# 			pdf(filename,8.27,11.69) # A4 size in inches
# 		} else {
# 			pdf(filename,ncol*2,ysize)
# 		}
# 	}
# 	k <- 1
# 	while(k<=objSize){
# 		par(mar=c(0,4,0,0))
# 		plot=plot(NA,xlim=c(0,24*ncol+txt*(ncol-1)),ylim=c(0.5,ysize),type="n",axes=FALSE,ann=FALSE)
# 		height <- ysize
# 		while(height>=1){
# 			if(k<=objSize){
# 				height <- height-max(minSpace,round(ord[k]+0.1,digits=1))
# 			}
# 			if(height>=1){
# 				#         cat("ligne ",count,"\n")
# 				count <- count+1
# 				for(j in 0:ncol-1){
# 					if(k<=objSize){
# 						#             cat(k,"\n")
# 						# Draw axis segments
# 						#             if(!colored){
# 						segments(0+j*(24+txt),height+1,24+j*(24+txt),height+1,lwd=0.5)
# 						#             }
# 						segments(0+j*(24+txt),height+1+0.05,0+j*(24+txt),height+1-0.05,lwd=0.5)
# 						segments(7+j*(24+txt),height+1+0.05,7+j*(24+txt),height+1-0.05,lwd=0.5)
# 						segments(14+j*(24+txt),height+1+0.05,14+j*(24+txt),height+1-0.05,lwd=0.5)
# 						segments(21+j*(24+txt),height+1+0.05,21+j*(24+txt),height+1-0.05,lwd=0.5)
# 						
# 						# Draw axis names
# 						# Version without protNames
# 						#             lab <- paste0(protAcc[k],"\n",geneName[k])
# 						# Version with protNames
#   					if(!is.null(names)){
#               pName <- protName[k]
#   						if(nchar(toString(pName))>30){
#   							pName <- paste0(substring(pName,1,15),"\n",substring(pName,16,min(27,nchar(toString(pName)))),"...")
#   						} else if(nchar(toString(pName))>15) {
#   							pName <- paste0(substring(pName,1,15),"\n",substring(pName,16,nchar(toString(pName))))
#   						}
#   						lab <- paste0(protAcc[k],"-",geneName[k],"\n",pName)
#   						
#   						axis(side=2, las=1, at=height+1, pos=j*(24+txt), labels=lab,cex.axis=0.4)
# 					  }
# 						
# 						# Draw plotnoise
# 						if (plotnoise) {
# 							.drawShading(obj[[k]]$targets$x + j*(24+txt), obj[[k]]$targets$pmean + height
# 										 , obj[[k]]$targets$pstd, color, obj[[k]]$targets$noisestd, uniformvar=uniformvar)
# 						} else {
# 							.drawShading(obj[[k]]$targets$x + j*(24+txt), obj[[k]]$targets$pmean + height
# 										 , obj[[k]]$targets$pstd, color, uniformvar=uniformvar)
# 						}
# 						# Draw plotmean
# 						if(colored){
# 							ymax = max(1.3,max(obj[[k]]$targets$npupper95))
# 							ymin = min(1/1.3,min(obj[[k]]$targets$nplower95))
# 							ylim = c(ymin,ymax)
# 							.drawRGB(obj[[k]]$targets$x + j*(24+txt), obj[[k]]$targets$pmean + height
# 										 , obj[[k]]$targets$pstd,lwd=lwd,ylimits=ylim,shift=height)
# 						}
# 						lines(obj[[k]]$targets$x + j*(24+txt), obj[[k]]$targets$pmean + height,
# 							  col=colorLine, lwd=lwd)
# 						k <- k+1
# 					} else {
# 						height <- 0
# 						k <- objSize+1
# 					}
# 					
# 				}
# 				if(k<=objSize){
# 					height <- height-max(minSpace,round(ord[k-ncol]+0.1,digits=1))
# 				}
# 			}
# 		}
# 		if(!file) { k <- objSize+1 }
# 	}
# 	if(file){
# 		dev.off()
# 	}
# }
# 
# 
# # Create a set of gppack plots. Either the whole data in a PDF (several pages) or just the first page in R.
# plot.gppackset <- function(obj, names=NULL, nrow=5, ncol=7, file=TRUE, filename=NULL, 
# 						   		plotdata=FALSE, plotnoise=FALSE) {
# 	N = length(obj)
# 	nx = round(0.8*sqrt(N))
# 	ny = ceiling(N/nx)
# 	
# 	if (!is.null(filename)) {
# 		pdf(filename, nx*1.5, ny*1.5)
# 	}
# 	get( getOption( "device" ) )()
# 	split.screen(figs = c(nx,ny))
# #	par(mfrow=c(nx,ny))
# 	for (i in 1:N) {
# 		screen(i)
# 		plot(obj[[i]])
# 	}
# 	close.screen(all=TRUE)
# 	if (!is.null(filename)) {
# 		dev.off()
# 	}
# 	
# 	k <- 1
# 	while(k<=N){
# 		split.screen(c(nx,ny))
# 		for(i in 0:(nx-1)){
# 			for(j in 1:ny){
# 				if(k<=N){
# 					screen(i*ny+j)
# 					plot(obj[[k]])
# 					if(length(names)==N){
# 						mtext(names[k],side=1)
# 					}
# 					k <- k+1
# 				}
# 			}
# 		}
# 		close.screen(all=T)
# 	}	
# }
# 
# 
# 
# ###############################################
# #################### QUERY ####################
# ###############################################
# 
# # Query the plots and informations of a gene
# query <- function(target, gdata, pdata, ggps, pgps, gGO=NULL, pGO=NULL,colored=F){
# 	grow = match(target, gdata[,"symbol"])
# 	prow = match(target, pdata[,"symbol"])
#   if(is.na(prow)&&is.na(grow)){
#     prow = match(target, pdata[,"uniprot"])
#   }
# 	
# 	hgncID <- geneName <- refseq <- protAccess <- protName <- ""
# 	if(!is.na(grow)){
# 		hgncID <- gdata[grow,"hgnc"]
# 		geneName <- gdata[grow,"symbol"]
# 		refseq <- gdata[grow,"refseq"]
# 	}
# 	if(!is.na(prow)){
# 		hgncID <- pdata[prow,"hgnc"]
# 		geneName <- pdata[prow,"symbol"]
# 		protAccess <- pdata[prow,"uniprot"]
# 		protName <- pdata[prow,"prot_name"]
# 	}
# 	
# 	if(geneName==""){
# 		cat("NO RESULT")
# 	}else{
# 		mainSplit <- rbind(c(0,1,0.7,1),c(0,1,0,0.7))
# 		split.screen(mainSplit)
# 		
# 		########## INFOS
# 		screen(1)
# 		par(mar=c(0,0,0,0))
# 		plot(NA,xlim=c(1,10),ylim=c(4,0.5),axes=F,ann=F)
# 		txt <- paste0("[GENE] ",geneName)
# 		text(1,0.8,txt,pos=4,cex=1.2)
# 		if(!is.na(grow)){
# 			txt <- paste0(hgncID,"  -  ",refseq)
# 		} else {
# 			txt <- hgncID
# 		}
# 		text(1,1.3,txt,pos=4)
# 		if(nchar(toString(protName))>40){
# 			#sub <- paste("-",substr(toString(protName),41,nchar(toString(protName))))
# 			#protName <- paste(substr(toString(protName),1,40),"-")
# 			protName <- paste0(substr(toString(protName),1,40),"-\n-"
# 							   ,substr(toString(protName),41,nchar(toString(protName))))
# 		} #else { sub <- "" }
# 		if(!is.na(prow)){
# 			txt <- paste0("[PROTEIN] ",protName)
# 		} else {
# 			txt <- ""
# 		}
# 		text(6,0.8,txt,pos=4,cex=1.2)
# 		text(6,1.3,protAccess,pos=4)
# 		
# 		if(!is.na(grow) && !is.null(gGO)){ # GENE GO
# 			txt <- "Process :"
# 			text(1,1.7,txt,pos=4,cex=0.85)
# 			txt <- .wrap.it(toString(gGO[[grow]][["BPterm"]]),170)
# 			if(nchar(txt)>480){ txt <- paste0(substring(txt,1,480),"...")}
# 			text(1,2.2,txt,pos=4,cex=0.85)
# 			txt <- "Function :"
# 			text(1,2.7,txt,pos=4,cex=0.85)
# 			txt <- .wrap.it(toString(gGO[[grow]][["MFterm"]]),170)
# 			if(nchar(txt)>300){ txt <- paste0(substring(txt,1,300),"...")}
# 			text(1,3.1,txt,pos=4,cex=0.85)
# 			txt <- "Component :"
# 			text(1,3.6,txt,pos=4,cex=0.85)
# 			txt <- .wrap.it(toString(gGO[[grow]][["CCterm"]]),170)
# 			if(nchar(txt)>300){ txt <- paste0(substring(txt,1,300),"...")}
# 			text(1,3.9,txt,pos=4,cex=0.85)
# 		} else if(!is.na(prow) && !is.null(pGO)){ # PROT GO
# 			txt <- "Process :"
# 			text(1,1.7,txt,pos=4,cex=0.85)
# 			txt <- .wrap.it(toString(pGO[[prow]][["BPterm"]]),170)
# 			if(nchar(txt)>480){ txt <- paste0(substring(txt,1,480),"...")}
# 			text(1,2.2,txt,pos=4,cex=0.85)
# 			txt <- "Function :"
# 			text(1,2.7,txt,pos=4,cex=0.85)
# 			txt <- .wrap.it(toString(pGO[[prow]][["MFterm"]]),170)
# 			if(nchar(txt)>300){ txt <- paste0(substring(txt,1,300),"...")}
# 			text(1,3.05,txt,pos=4,cex=0.85)
# 			txt <- "Component :"
# 			text(1,3.6,txt,pos=4,cex=0.85)
# 			txt <- .wrap.it(toString(pGO[[prow]][["CCterm"]]),170)
# 			if(nchar(txt)>300){ txt <- paste0(substring(txt,1,300),"...")}
# 			text(1,3.95,txt,pos=4,cex=0.85)
# 		}
# 		
# 		########## PLOTS 
# 		screen(2)
# 		if(!is.na(grow) && !is.na(prow)){
# 			plotscr <- split.screen(c(1,2))
# 			screen(plotscr[1])
# 			plot(ggps[[grow]])
# 			screen(plotscr[2])
# 			plot(pgps[[prow]],colored=colored)
# 		} else if(!is.na(grow)) {
# 			plot(ggps[[grow]])
# 		} else if(!is.na(prow)) {
# 			plot(pgps[[prow]],colored=colored)
# 		}
# 		close.screen(all=T)
# 	}
# }
# 
# 
# 
# 
# plotemll = function(res) {
# 	plot(timepoints, rep(1,175), ylim=c(-3,1.5), t='l')
# 	lines(timepoints, res$nullmodel$targets$emll, col='green', lty=2)
# 	lines(timepoints, res$casemodel$targets$emll, col='blue', lty=2)
# 	lines(timepoints, res$ctrlmodel$targets$emll, col='red', lty=2)
# 	lines(timepoints, 0.5*res$ctrlmodel$targets$emll+0.5*res$casemodel$targets$emll-res$nullmodel$targets$emll, col='orange', lty=2)
# 	
# 	ctrl.f.mean = res$ctrlmodel$targets$pmean - mean(res$ctrlmodel$y)
# 	ctrl.k.m = res$ctrlmodel$ekernel
# 	ctrl.f.cov = res$ctrlmodel$cov
# 	ctrl.k.y = res$ctrlmodel$kernel
# 	ctrl.noise = res$ctrlmodel$target$noisestd
# 	ctrl.emll <- -0.5*ctrl.f.mean^2/diag(ctrl.k.m) -0.5*log(diag(ctrl.k.m)) -0.5*log(2*pi)
# 	
# 	case.f.mean = res$casemodel$targets$pmean - mean(res$casemodel$y)
# 	case.k.m = res$casemodel$ekernel
# 	case.f.cov = res$casemodel$cov	
# 	case.emll <- -0.5*case.f.mean^2/diag(case.k.m) -0.5*log(diag(case.k.m)) -0.5*log(2*pi)
# 	null.f.mean = res$nullmodel$targets$pmean - mean(res$nullmodel$y)
# 	null.k.m = res$nullmodel$ekernel
# 	null.f.cov = res$nullmodel$cov
# 	null.emll <- -0.5*null.f.mean^2/diag(null.k.m) -0.5*log(diag(null.k.m)) -0.5*log(2*pi)
# 	
# 	
# 	lines(timepoints, ctrl.f.mean, lty=2)
# 	lines(timepoints, diag(ctrl.k.m))
# 	
# 	plot(timepoints, rep(1,175), ylim=c(-2,1), t='l')
# 	lines(timepoints, null.emll, col='green')
# 	lines(timepoints, case.emll, col='blue')
# 	lines(timepoints, ctrl.emll, col='red')
# 	lines(timepoints, 0.5*case.emll+0.5*ctrl.emll-null.emll, col='orange')
# }	
# 


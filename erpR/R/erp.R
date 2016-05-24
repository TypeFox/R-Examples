####
# erp ##

erp <- function(el, startmsec=-200, endmsec = 1200, smo = NULL, interval = c(startmsec, endmsec), # here parameters for plot
xlim = c(startmsec, endmsec), ylim = c(-6,6), col = "black", lwd = 1, xlab = "", ylab = "", main = "", frame.plot = TRUE,
draw.xaxis = TRUE, draw.yaxis = TRUE, type = "l", x0 = TRUE,
x.tick = seq(-200, 1200, 200), x.labels=x.tick, x.pos = NA, #here parameters for x.axis
x.outer = F, x.font = NA, x.lty = "solid", x.lwd = 1, x.lwd.ticks = 1, x.col = NULL, 
x.col.ticks = NULL, x.hadj = NA, x.padj = NA, x.tcl = -0.5, x.tick.both= FALSE, x.cex = 1,
y.tick = seq(-6,6,2), y.labels = y.tick,  y.pos = NA,   # here parameters for y.axis
y.outer = F, y.font = NA, y.lty = "solid", y.lwd = 1, y.lwd.ticks = 1, y.col = NULL, 
y.col.ticks = NULL, y.hadj = NA, y.padj = NA, y.tcl = -0.5, y.tick.both= FALSE, y.cex = 1, ...)
 
# NOTICE!!!! in the parameters on y-axis compared to the erp.yaxis function, the only difference is in the default of y.pos, that is startmsec 
# (NA would work bad, because another axis is drawn before)

# the three dots refer to plot
{
	
	# first I calculate the msectopoints on the WHOLE length of the data
	lengthwhole=length(el)	
	
	startpoint=msectopoints(xlim[1], lengthwhole, startmsec, endmsec)
	endpoint=msectopoints(xlim[2], lengthwhole, startmsec, endmsec)
	
	# determine the waveform ot plot according to interval
	plotinterval=msectopoints(interval[1], lengthwhole, startmsec, endmsec): msectopoints(interval[2], lengthwhole, startmsec, endmsec )
	
	
	if (!is.null(smo)){
		el=smooth.spline(el, spar = smo)$y
	}
	
	### PLOT WAVEFORM
	# notice the default xlim is from 1 to the whole length of the data (not in intervals)
	
	plot(plotinterval, el[plotinterval], xlim =c(startpoint, endpoint), ylim = ylim, col = col, lwd = lwd, xlab = xlab, ylab = ylab, main = main, frame.plot = frame.plot, 
	axes=F, type = type, ...)
	
	# the only default plotting is a segment indicating the 0.
	if (x0 == TRUE){
		x0.l=as.numeric(ylim)/10
		zeropoint=msectopoints(0, lengthwhole, startmsec, endmsec)
		segments(zeropoint, x0.l, zeropoint, -x0.l, lwd=1.5)
	}
	
	
	# PLOT X.AXIS calling erp.xaxis function
	# calculate default for axis according to xlim provided
	
	
	if (draw.xaxis==TRUE){
	erp.xaxis(lengthwhole, startmsec = startmsec, endmsec = endmsec, x.labels = x.labels, x.pos=x.pos, x.tick=x.tick,
	 x.outer=x.outer, x.font=x.font, x.lty=x.lty, x.lwd=x.lwd, x.lwd.ticks=x.lwd.ticks, x.col=x.col, x.col.ticks=x.col.ticks, x.hadj=x.hadj, x.padj=x.padj, 
	 x.tcl, x.tick.both=x.tick.both, x.cex=x.cex)
	 }
	
	 # PLOT y.axis calling erp.yaxis function
	 # calculate default for axis according to ylim provided
	 

	if (draw.yaxis==TRUE){
	erp.yaxis(lengthwhole, startmsec = startmsec, endmsec = endmsec, y.labels = y.labels, y.pos=y.pos, y.tick=y.tick,
	 y.outer=y.outer, y.font=y.font, y.lty=y.lty, y.lwd=y.lwd, y.lwd.ticks=y.lwd.ticks, y.col=y.col, y.col.ticks=y.col.ticks, 
	 y.hadj=y.hadj, y.padj=y.padj, y.tcl, y.tick.both=y.tick.both, y.cex=y.cex)
	}
	
}
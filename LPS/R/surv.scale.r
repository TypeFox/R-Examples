# Draws a survival scale for legend
# Author : Sylvain Mareschal <maressyl@gmail.com>
surv.scale <- function(time, event, eventColors=c("#000000", "#CCCCCC"), censColors=c("#FFFFEE", "#FFDD00")) {
	# Layout
	layout(matrix(1:2, ncol=1))
	par(oma=c(1,0,0,0))
	
	# With events
	par(mar=c(2.5,1,1,1))
	plot(x=NA, y=NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
	rect(xleft=(0:255)/256, xright=(1:256)/256, ybottom=0, ytop=1, border=NA, col=rgb(colorRamp(eventColors)((0:255)/255), maxColorValue=255))
	axis(side=1, at=0:1, labels=c(0, max(time[ event ], na.rm=TRUE)))
	mtext(side=1, text="until metastasis", line=1)
	box()
	
	# Censored
	par(mar=c(3,1,0.5,1))
	plot(x=NA, y=NA, xlim=0:1, ylim=0:1, xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i")
	rect(xleft=(0:255)/256, xright=(1:256)/256, ybottom=0, ytop=1, border=NA, col=rgb(colorRamp(censColors)((0:255)/255), maxColorValue=255))
	axis(side=1, at=0:1, labels=c(0, max(time[ !event ], na.rm=TRUE)))
	mtext(side=1, text="until censoring", line=1)
	box()
	
	mtext(side=1, text="Follow up (months)", outer=TRUE)
}


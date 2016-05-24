## ==============================
## PLot of the state distribution
## ==============================

plot.stslist.statd <- function(x, type="d", cpal=NULL, ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, xtstep=NULL, 
	cex.plot=1, space=0, ...) {

	n <- attr(x,"nbseq")
	weighted <- attr(x, "weighted")
	if (weighted) {wlab <- "weighted "}
	else {wlab <- NULL}

	if (is.null(xtlab)) {xtlab <- attr(x,"xtlab")}

	if (is.null(xtstep)) {
		if (!is.null(attr(x,"xtstep"))) {xtstep <- attr(x,"xtstep")} 
		## For sequence objects created with previous versions
		else {xtstep <- 1}
	}

	seql <- length(xtlab)

	## State distribution plot	
	if (type=="d") {
		if (is.null(cpal))
			cpal <- attr(x,"cpal")

		if (is.null(ylab)) 
			ylab <- paste("Freq. (",wlab,"n=",round(n,2),")",sep="")

		barplot(x$Frequencies,
			space=space,
			axes=FALSE,
			axisnames=FALSE,
			## cex.axis=cex.plot,
			## cex.names=cex.plot,
			col=cpal,
			ylab=ylab,
			...)

		## Plotting the x axis
		x.lab.pos <- NULL
		tpos <- seq(1,seql, xtstep)

		for (p in tpos) {
			x.lab.pos <- c(x.lab.pos, (p-1)+((p-1)*space)+(0.5+space))
		}

		if (xaxis) {
			axis(1, at=x.lab.pos, labels=xtlab[tpos], pos=-0.02, cex.axis=cex.plot)
		}
	}
	## Entropy index plot
	else if (type=="Ht") {
		if (is.null(ylab)) 
			ylab <- paste("Entropy index (",wlab,"n=",round(n,2),")",sep="")
	
		plot(x$Entropy, 
			col="blue",
			## frame.plot=TRUE,
			type="l",
			lwd=3.5, 
			lty="solid", 
			axes=FALSE,
			ylim=0:1,
			ylab=ylab,
			xlab=NA,
			...)

		## Plotting the x axis
		if (xaxis) {
			tpos <- seq(1,seql, xtstep)
			axis(1, at=tpos, labels=xtlab[tpos], pos=-0.02, cex.axis=cex.plot)
		}
	}

	##
	if (is.null(yaxis) || yaxis)	
		axis(2, cex.axis=cex.plot)
}

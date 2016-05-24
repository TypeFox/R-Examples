## ================================
## PLot of the sequences frequency
## ================================

plot.stslist.freq <- function(x, cpal=NULL, missing.color=NULL, pbarw=TRUE, 
	ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, xtstep=NULL, cex.plot=1, ...) {

	n <- attr(x,"nbseq")
	weighted <- attr(x, "weighted")
	if (weighted) {wlab <- "weighted "}
	else {wlab <- NULL}
	statl <- attr(x,"alphabet")
	nr <- attr(x,"nr")

	if (is.null(xtlab))
		xtlab <- attr(x,"names")

	if (is.null(xtstep)) {
		if (!is.null(attr(x,"xtstep"))) {xtstep <- attr(x,"xtstep")} 
		## For sequence objects created with previous versions
		else {xtstep <- 1}
	}

	seql <- length(xtlab)

	if (is.null(cpal))
		cpal <- attr(x,"cpal")

	## Checking for missing states
	if (any(x==nr)) {
		if (is.null(missing.color)) missing.color <- attr(x,"missing.color")
		cpal <- c(cpal, missing.color)
		statl <- c(statl, nr)
	}

	if (is.null(ylab)) {
		if (yaxis==TRUE || yaxis=="cum")
			ylab <- paste("Cum. % freq. (",wlab,"n=",round(n,2),")",sep="")
		else if (yaxis=="pct")
			ylab <- paste("% freq. (",wlab,"n=",n,")",sep="")		
	}

	## Storing the optional parameters in a list
	olist <- list(...)

	tlim <- nrow(x)

	## 
	seqbar <- apply(x,1, seqgbar, seql=seql, statl=statl)

	table <- attr(x,"freq")

	if (pbarw==TRUE) barw=table$Percent 
	else barw=1

	## The plot
	barplot(seqbar,col=cpal, width=barw,
		ylab=ylab,
		horiz=TRUE,
		axes=FALSE,
		axisnames=FALSE,
		...)
	
	## Plotting the x axis
	if (xaxis) {
		tpos <- seq(1, seql, xtstep)
		axis(1, at=tpos-0.5, labels=xtlab[tpos], cex.axis=cex.plot)
	}

	## Plotting the y axis
	if ("space" %in% names(olist)) space <- olist[["space"]]
	else space <- 0.2

	if (yaxis==TRUE || is.null(yaxis) || yaxis=="cum") {
		y.lab <- paste(c(0, round(sum(table$Percent),1)),"%",sep="")
		y.tick <- TRUE

		if (!pbarw)
			y.lab.pos <- c(space,(tlim+(tlim*space)))
		else
			y.lab.pos <- c(space*mean(barw),(sum(barw)+(tlim*space*mean(barw))))
	}
	## Percentage frequency of each sequence 
	else if (yaxis=="pct") {
		y.lab <- round(table$Percent,1)
		y.tick <- FALSE

		if (!pbarw) {
			y.lab.pos <- 0.7
			for (p in 2:length(y.lab))
				y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*space)+0.7)
			} 
		else { 
			y.lab <- y.lab[table$Percent>=0.5]

			y.lab.pos <- (table$Percent[1]/2)+1
			sep <- space*mean(table$Percent)

			for (p in 2:length(y.lab))
				y.lab.pos <- c(y.lab.pos, sum(table$Percent[1:p])+(p*sep)-table$Percent[p]/2)
			} 
	}

	if (yaxis==TRUE || yaxis=="cum" || yaxis=="pct")
		axis(2, at=y.lab.pos, 
			labels=y.lab, 
			tick=y.tick,
			## mgp=c(1.5,1,0), 
			las=1, 
			cex.axis=cex.plot)
}
	

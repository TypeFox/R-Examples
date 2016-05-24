## =============================
## Plot the modal state sequence
## =============================

plot.stslist.modst <- function(x, cpal=NULL,
	ylab=NULL, yaxis=TRUE, xaxis=TRUE, xtlab=NULL, xtstep=NULL, cex.plot=1, ...) {

	seql <- ncol(x)
	statl <- attr(x,"alphabet")
	n <- attr(x, "nbseq")
	nr <- attr(x,"nr")

	if (is.null(cpal)) {cpal <- attr(x,"cpal")}

	## Adding an entry for missing in the legend
	if (any(x==nr)) {
		missing.color <- attr(x,"missing.color")
		cpal <- c(cpal, missing.color)
		statl <- c(statl, nr)
	}
	nbstat <- length(statl)

	weighted <- attr(x, "weighted")
	if (weighted) {wlab <- "weighted "}
	else {wlab <- NULL}

	if (is.null(xtlab)) {xtlab <- colnames(x)}
	if (is.null(xtstep)) {
		if (!is.null(attr(x,"xtstep"))) {xtstep <- attr(x,"xtstep")}
		## For sequence objects created with previous versions
		else {xtstep <- 1}
	}

	if (is.null(ylab)) ylab <- paste("State freq. (",wlab,"n=",round(n,2),")",sep="")

	## ============================
	## Max distance for axis limits
	## ============================
	mod <- attr(x, "Frequencies")
	prof.freq <- matrix(0, nrow=nbstat, ncol=seql)

	## Preparing the matrix for plot
	for (i in 1:seql) {
		smax <- which(statl==x[,i])
		prof.freq[smax,i] <- mod[i]
	}

	## Frequency of the representative sequence
	nbrep <- attr(x,"Occurrences")
	ctfreq <- round((nbrep/n)*100,1)
	txt <- paste("Modal state sequence (",nbrep," occurrences, freq=", ctfreq ,"%)", sep="")

	barplot(prof.freq,
		space=0,
		## mgp=c(2.5,0.6,0),
		cex.names=cex.plot,
		ylim=c(0,1.2),
		col=cpal,
		## main=title,
		axisnames=FALSE,
		ylab=ylab,
		axes=FALSE,
		...)

	text(seql/2, 1.1, txt,
		cex=cex.plot)

	## Plotting the x axis
	if (xaxis) {
		tpos <- seq(1,seql, xtstep)
		axis(1, at=tpos-0.5, labels=xtlab[tpos], pos=-0.02,
		# mgp=c(3,0.5,0),
		cex.axis=cex.plot)
	}

	## Axis for the state frequencies
	if (yaxis)
		axis(2, at=seq(0,1.0,0.25), labels=c("0","0.25",".5","0.75","1"),
			las=2, cex.axis=cex.plot)

}

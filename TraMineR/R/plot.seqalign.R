## =============================
## PLot of STS sequence objects
## =============================

plot.seqalign <- function(x,
	cpal=NULL, missing.color=NULL, ylab=NULL, yaxis=TRUE, xaxis=TRUE, ytlab=NULL, ylas=0, xtlab=NULL, cex.plot=1, ...) {

	showop <- "bars"
	n <- 2
	seql <- length(x$seq1)
	statl <- alphabet(x$stsseq)
	nr <- attr(x$stsseq,"nr")

	op <- x$operation
	cost <- x$cost

	##
	if (is.null(cpal)) {cpal <- attr(x$stsseq,"cpal")}

	## Adding an entry for missing in the legend
	if (any(x$stsseq==nr)) {
		if (is.null(missing.color)) missing.color <- attr(x$stsseq,"missing.color")
		cpal <- c(cpal, missing.color)
		statl <- c(statl, nr)
	}

	cpal <- c(cpal, "white")
	statl <- c(statl, "-")
	

	## Storing the optional parameters in a list
	olist <- list(...)

	ssamp <- rbind(x$seq2, x$seq1)
	seqbar <- apply(ssamp, 1, seqgbar, statl=statl, seql=seql)

	seqop <- apply(matrix(op,nrow=1), 1, seqgbar, statl=c("D","E","I","S"), seql=seql)
	dummy <- rep(0,length(seqop))

	seqins <- cbind(dummy, dummy, seqop)
	seqsub <- cbind(dummy, dummy, dummy, seqop)
	seqeq <- cbind(dummy, dummy, dummy, dummy, seqop)

	ylab <- "Alignment"

	## The PLot
	barplot(seqbar,col=cpal, ylim=c(0,4),
		## ylab=ylab,
		horiz=TRUE,
		yaxt="n",
		axes=FALSE,
		las=1,
		...
	)

	text(seql/2, 2.5, paste("Alignment (cost:", sum(x$cost),")"))

	albgd <- "grey90"

	barplot(seqins, add=TRUE, col=c("grey50",albgd,"grey50",albgd), width=c(1,1,.3),
		## horiz=TRUE,
		yaxt="n",
		axes=FALSE, horiz=TRUE, space=.35, ...
		)

	barplot(seqsub, add=TRUE, col=c(albgd,albgd,albgd,"red"), width=c(1,1,1.1,.4),
		## horiz=TRUE,
		yaxt="n",
		axes=FALSE, horiz=TRUE, space=0, ...
		)


	barplot(seqeq, add=TRUE, col=c(albgd,"green",albgd,albgd), width=c(1,1,1.1,.3,.3),
		## horiz=TRUE,
		yaxt="n",
		axes=FALSE, horiz=TRUE, space=0, ...
		)

	rect(0, 2.8, seql, 3.7)

	text(seql/2, 3.9,
		paste("Operations (", sum(op=="S"), " subst., ", sum(op %in% c("I","D"))," indels)", sep=""))

	if (showop=="symbol") {
		for (i in 1:length(op)) {
			if (op[i]=="I") {
				arrows(i-0.55, 1.3, i-0.45, 1.3, length=0.1, angle=30, col="red", lwd=2)
				text(i-0.5,1.25, cost[i], col="red", cex=0.9, font=2)
			}
			else if (op[i]=="D") {
				arrows(i-0.8, 1.3, i-0.7, 1.3, code=1, length=0.1, angle=30, col="red", lwd=2)
				text(i-0.5,1.25, cost[i], col="red", cex=0.9, font=2)
			}
			else if (op[i]=="S") {
				arrows(i-0.5, 1.29, i-0.5, 1.31, code=1, length=0.1, angle=30, col="blue", lwd=2)
				text(i-0.5,1.25, cost[i], col="blue", cex=0.9, font=2)
			}
		}
	}

	## Plotting the x axis
	if (xaxis) {
	axis(1, at=1:seql-0.5, labels=1:seql,
		## mgp=c(3,0.5,0),
		cex.axis=cex.plot)
	}


	## Plotting the y axis
	if (is.null(yaxis) || yaxis) {
		if ("space" %in% names(olist)) sp <- olist[["space"]]
		else sp <- 0.2
	
		y.lab.pos <- sp+0.5
		y.lab.pos <- c(y.lab.pos, 1+(1*sp)+(0.5+sp))

		if (is.null(ytlab)) {ytlab <- paste("seq",2:1, sep="")}
		## else if (ytlab=="id") {ytlab <- rownames(x)[tlim]}

		axis(2, at=y.lab.pos, mgp=c(1.5,0.5,0), labels=ytlab, las=ylas, tick=FALSE, cex.axis=cex.plot)

		lab.op.pos <- c(2.7+sp, 3.0+sp, 3.3+sp)
		axis(2, at=lab.op.pos, mgp=c(1.5,0.5,0), labels=c("IND","SUB","EQU"),
			las=2, tick=FALSE, cex.axis=cex.plot)
		
	}

}


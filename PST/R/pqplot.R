## Probability distribution of a node and its parents
## and outcomes of the gain function 

setMethod("pqplot", signature=c(object="PSTf", data="stslist"), 
	def <- function(object, data, cdata, L, stcol, plotseq=FALSE, ptype="b", 
	cex.plot=1, space=0, measure="prob", pqmax, seqscale, ...) {

	oolist <- list(...)

	if (missing(L)) { L <- length(object)-1 }
	
	message(" [>] computing prob., max. length=", L)
	
	if (!missing(cdata)) {
		prob <- suppressMessages(predict(object, data=data, cdata=cdata, L=L, decomp=TRUE))
	} else {
		prob <- suppressMessages(predict(object, data=data, L=L, decomp=TRUE))
	}

	## Number of predicted symbols (used instead of sequence length)
	sl <- rowSums(!is.na(prob))

	if (measure=="logloss") {
		prob <- -log(prob, base=2)
		if (missing(pqmax)) { pqmax <- ceiling(max(prob)) }
		pmean <- sum(prob)/sl
		ytstep <- 1
		ylab <- "Log-loss"
	} else {
		if (missing(pqmax)) { pqmax <- 1 }
		pmean <- exp(log(rowProds(prob))/sl)
		ytstep <- 0.2
		ylab <- "Prob"
	}

	poff <- 0

	if (plotseq) {
		if (missing(cdata)) { cdata <- data }

		c.A <- alphabet(cdata)
		c.cpal <- cpal(cdata)

		if (any(cdata==attr(cdata, "nr"))) {
			c.A <- c(c.A, attr(cdata, "nr"))
			c.cpal <- c(c.cpal, attr(cdata, "missing.color"))
		}

		tmp <- TraMineR:::seqgbar(as.matrix(cdata), seql=sl, statl=c.A)

		## Plotting path
		barw <- 1
		if (missing(seqscale)) { seqscale <- pqmax/3 }
		seqpsep <- 0.1*pqmax
		poff <- seqscale+seqpsep

		tmp <- matrix(tmp, nrow=length(c.A))
		tmp <- tmp*seqscale

		barplot(tmp, col=c.cpal, width=barw,
			## ylab=ylab,
			## xlim=c(0,(sl+1)),
			ylim=c(0,(poff+pqmax)),
			horiz=FALSE,
			axes=FALSE,
			axisnames=FALSE,
			space=space,
			...)
		if (missing(stcol)) { stcol="grey" }

	}

	## Plotting probability distributions
	if (ptype=="b") {
		if (missing(stcol)) {
			A <- alphabet(data)
			cpal <- cpal(data)
			nr <- attr(data, "nr")

			if (any(data==attr(data, "nr"))) {
				A <- c(A, nr)
				cpal <- c(cpal, attr(data, "missing.color"))
			}

			tmp <- TraMineR:::seqgbar(as.matrix(data), seql=sl, statl=A)
			tmp[tmp==1] <- tmp[tmp==1]*prob
			tmp <- matrix(tmp, nrow=length(A))
			stcol <- cpal
		} else {
			tmp <- prob
		}

		barplot(tmp, col=stcol, offset=poff, add=plotseq, ylim=c(0,(poff+pqmax)), 
			space=space, axes=FALSE, axisnames=FALSE, ylab=NULL, ...)
		abline(h=pmean+poff, col="red")

	} else if (ptype=="l") {
		if (missing(stcol)) { stcol <- "blue" }

		xt <- if (ptype=="s") { 0:(length(prob)-1) } else { 0.5:(length(prob)-0.5) }

		if (plotseq) {
			lines(xt, prob+poff, col=stcol, type=ptype)
		} else {
			plot(xt, prob+poff, col=stcol, type=ptype)
		}
	}

	tpos <- seq(1, sl, by=attr(data, "xtstep"))
	tlab <- colnames(data)[tpos]
	tpos <- tpos + (tpos*space)

	axis(1, at=tpos-0.5, labels=tlab, pos=-0.04)

	plabpos <- seq(from=poff, to=(poff+pqmax), by=ytstep)
	plab <- plabpos-poff

	axis(2, at=plabpos, 
		labels=plab, 
		## las=2, 
		cex.axis=cex.plot)

	mtext(ylab, side=2, at=(min(plabpos)+max(plabpos))/2, line=3)

}
)

## =============================
## PLot of STS sequence objects
## =============================

plot.stslist <- function(x, tlim=NULL, weighted=TRUE, sortv=NULL,
	cpal=NULL, missing.color=NULL, ylab=NULL, yaxis=TRUE, xaxis=TRUE, ytlab=NULL, ylas=0, xtlab=NULL, xtstep=NULL, cex.plot=1, ...) {

	n <- nrow(x)
	seql <- ncol(x)
	statl <- alphabet(x)
	nr <- attr(x,"nr")

	if (is.null(xtlab))
		xtlab <- colnames(x)

	if (is.null(xtstep)) {
		if (!is.null(attr(x,"xtstep"))) {xtstep <- attr(x,"xtstep")}
		## For sequence objects created with previous versions
		else {xtstep <- 1}
	}

	## Range
	if (is.null(tlim)) {
		if (n>=10) tlim <- 1:10
		else tlim=1:n
	}
	else if (tlim[1]==0)
			tlim <- 1:n
	else if (max(tlim) > n)
			tlim <- 1:n
	
	## Sorting
	if (!is.null(sortv)) {
		if (length(sortv)==1 && sortv %in% c("from.start", "from.end")) {
        		end <- if (sortv=="from.end") { max(seqlength(x)) } else { 1 }
        		beg <- if (sortv=="from.end") { 1 } else { max(seqlength(x)) }

			sortv <- do.call(order, as.data.frame(x)[,end:beg])
			x <- x[sortv,]
		} else if (length(sortv)!=n) {
			stop(call.=FALSE, "sortv must contain one value for each row in the sequence object")
		} else {
			if (is.factor(sortv)) { sortv <- as.integer(sortv) }
			x <- x[order(sortv),]
		}

		sortlab <- paste(", sorted")
		
	} else { sortlab <- NULL }

	##
	if (is.null(cpal))
		cpal <- attr(x,"cpal")

	## Adding an entry for missing in the legend
	if (any(x==nr)) {
		if (is.null(missing.color)) missing.color <- attr(x,"missing.color")
		cpal <- c(cpal, missing.color)
		statl <- c(statl, nr)
	}

	## Storing the optional parameters in a list
	olist <- list(...)

	ssamp <- x[tlim,]
	seqbar <- apply(ssamp, 1, seqgbar, statl=statl, seql=seql)

	## WEIGHTS
	## Weights
	weights <- attr(x, "weights")

	if (!weighted || is.null(weights)) {
		weights <- rep(1.0, nrow(x))
	}
	## Also takes into account that in unweighted sequence objects created with
	## older TraMineR versions the weights attribute is a vector of 1
	## instead of NULL
	if (all(weights==1))
		weighted <- FALSE

	if (weighted) {wlab <- "weighted "}
	else {wlab <- NULL}

	if (is.null(ylab)) {
		ylab <- paste(length(tlim)," seq. ", "(", wlab,"n=", round(sum(weights),2),")",
			sortlab, sep="")
	}

	## The PLot
	barplot(seqbar,col=cpal, width=weights,
		ylab=ylab,
		horiz=TRUE,
		yaxt="n",
		axes=FALSE,
		las=1,
		...
	)

	## Plotting the x axis
	if (xaxis) {
		tpos <- seq(1,seql, xtstep)
		axis(1, at=tpos-0.5, labels=xtlab[tpos],
		# mgp=c(3,0.5,0),
		cex.axis=cex.plot)
	}

	## Plotting the y axis
	if (is.null(yaxis) || yaxis) {
		if ("space" %in% names(olist)) sp <- olist[["space"]]
		else sp <- 0.2

		idxmax <- length(tlim)
		
		if (!weighted) {
			y.lab.pos <- sp+0.5

			if (idxmax>1) {
				for (p in 2:idxmax) {
					y.lab.pos <- c(y.lab.pos, (p-1)+((p-1)*sp)+(0.5+sp))
				}
			}
		}
		else {
			y.lab.pos <- (weights[1]/2)+sp
			sep <- sp*mean(weights)

			if (idxmax>1) {
				for (p in 2:idxmax)
					y.lab.pos <- c(y.lab.pos, sum(weights[1:p])+(p*sep)-weights[p]/2)
			}
		}

		if (is.null(ytlab)) {ytlab <- tlim}
		else if (ytlab=="id") {ytlab <- rownames(x)[tlim]}

		axis(2, at=y.lab.pos, mgp=c(1.5,0.5,0), labels=ytlab, las=ylas, tick=FALSE, cex.axis=cex.plot)
	}

}


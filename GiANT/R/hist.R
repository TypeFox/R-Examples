hist.gsaResult <- function(x,
		signLevel = x$signLevel,
		subset = NULL,
		ask = FALSE,
		addLegend = TRUE,
		...){
	#print("GSA histogram")
	oldAsk <- par("ask")

	args <- list(...)

	if(is.null(subset)){
		subset <- 1:length(x$res.all)
	}

	##########################################
	# check graphical parameters for use with
	# list of results in x
	##########################################
	# main
	if(is.null(args$main)){
		mains <- names(x$res.all)
	}else if((length(args$main) != 1) && (length(args$main) != length(x$res.all))){
		stop("Incorrect length of main (length = ", length(args$main),
			").\nMust be 1 or equal to the number of analyzed gene sets or NULL (for default).", sep="")
	}else if(length(args$main) == 1){
		mains <- rep(args$main, length(x$res.all))
	}else{
		mains <- args$main
	}
	# xlab
	if(is.null(args$xlab)){
		xlabs <- rep("", length(x$res.all))
	}else if((length(args$xlab) != 1) && (length(args$xlab) != length(x$res.all))){
		stop("Incorrect length of xlab (length = ", length(args$xlab),
			").\nMust be 1 or equal to the number of analyzed gene sets or NULL (for no xlab).", sep="")
	}else if(length(args$xlab) == 1){
		xlabs <- rep(args$xlab, length(x$res.all))
	}else{
		xlabs <- args$xlab
	}
	# ylab
	if(is.null(args$ylab)){
		ylabs <- rep("", length(x$res.all))
	}else if((length(args$ylab) != 1) && (length(args$ylab) != length(x$res.all))){
		stop("Incorrect length of ylab (length = ", length(args$ylab),
			").\nMust be 1 or equal to the number of analyzed gene sets or NULL (for no ylab).", sep="")
	}else if(length(args$ylab) == 1){
		ylabs <- rep(args$ylab, length(x$res.all))
	}else{
		ylabs <- args$ylab
	}

	xlim <- NULL
	ylim <- NULL

	if(!is.null(args$xlim)){
		xlim <- args$xlim
	}
	if(!is.null(args$ylim)){
		ylim <- args$ylim
	}

	args <- args[-which(names(args) %in% c("main", "xlab", "ylab", "xlim", "ylim"))]
	##########################################
	for(i in subset){

		res.all <- x$res.all[[i]]
		pValue <- x$adjustedPValues[i]

		numSamples <- length(res.all$significanceValues$gssValues)

		if(is.null(xlim)){
			dat.min <- min(c(res.all$geneSetValues$gss, res.all$significanceValues$gssValues))
			dat.max <- max(c(res.all$geneSetValues$gss, res.all$significanceValues$gssValues))
			xlim <- c(dat.min,dat.max)
		}
		if(is.null(ylim)){
			ylim <- c(0,numSamples/6)
		}
		##########################################
		## plot histogram enrichment analyses
		##########################################
		par(ask=ask && i > 1 && dev.interactive())
		
		histObjects <- do.call(hist, args=c(
			list(
				x = c(res.all$geneSetValues$gss, res.all$significanceValues$gssValues),
				breaks = 50, 
				col="lightblue",
				main=paste(mains[i], " (p = ",round(pValue, digits=4),")",sep=""),
				xlab=xlabs[i],
				ylab=ylabs[i],
				xlim=xlim,
				ylim=ylim),
			args))

		# mark 5% and 95% quantile
		if(x$analysis$testAlternative == "greater"){
			quant <- quantile(res.all$significanceValues$gssValues,1-signLevel)
			th <- 1-signLevel
		}else if(x$analysis$testAlternative == "less"){
			quant <- quantile(res.all$significanceValues$gssValues,signLevel)
			th <- signLevel
		}

		abline(v=quant,
			col="blue",
			lwd=2,
			lty=3)
		# mark value of tested geneset
		abline(v=res.all$geneSetValues$gss,
			col="red",
			lwd=2,
			lty=1)
			
		# legend
		if(addLegend){
			tmp_t <- substitute(expression(t[tmp], "t*"), env = list(tmp = th))
			legend("topleft",
				legend=eval(tmp_t),
				lty=c(3,1),
				#bty="n",
				bg = "white",
				col=c("blue","red"))
		}
	}
	par(ask=oldAsk)

	return(invisible(histObjects))
}

plot.uncertaintyResult <- function(x,
		signLevel = x$signLevel,
		addLegend = TRUE,
		addMinimalStability = FALSE,
		...){

	args <- list(...)
	quant <- c(signLevel, 0.5, 1-signLevel)

	k <- sapply(x$uncertaintyEvaluations, "[[", 4)

	if(signLevel == x$signLevel){
		quantiles <- x$confidenceValues
	}else{
		quantiles <- t(sapply(x$uncertaintyEvaluations, function(i){
				quantile(i$gssValues, probs = quant)
			}))
		rownames(quantiles) <- rownames(x$confidenceValues)
	}

	##########################################
	# check graphical parameters for use with
	# list of results in x
	##########################################
	# main
	if(is.null(args$main)){
		myMain <- ""
	}else{
		myMain <- paste(args$main, "\n", sep ="")
	}
	# xlab
	if(is.null(args$xlab)){
		xlab <- "original geneSet genes (k)"
	}else{
		xlab <- args$xlab
	}
	# ylab
	if(is.null(args$ylab)){
		ylab <- "test statistic (t)"
	}else{
		ylab <- args$ylab
	}

	nullDistr <- quantile(x$originalGeneSetValues$res.all[[1]]$significanceValues$gssValues, probs = quant)
	t0 <- x$originalGeneSetValues$res.all[[1]]$geneSetValues$gss

	quantiles <- rbind(nullDistr, quantiles, rep(t0, length(quant)))

	newk <- c(0,k,1)

	if(addMinimalStability){
		a <- (x$originalGeneSetValues$res.all[[1]]$geneSetValues$gss - nullDistr[2])
		linEstimate <- (nullDistr[3]-nullDistr[1])/a
		myMain <- paste(myMain, "resampling stability:", round(x$uncertainty, digits = 2), "\n",
			"minimal stability:", round(linEstimate, digits = 2), sep = " ")
	}else{
		myMain <- paste(myMain, "estimated uncertainty:", x$uncertainty, sep = " ")
	}

	plot(y=quantiles[,1],
		x = newk,
		ylab = ylab,
		xlab = xlab,
		type = "o",
		main = myMain,
		col = "#045a8d",
		lwd = 2,
		ylim = range(quantiles)*c(0.975,1.025))
	lines(quantiles[,2], x = newk, type = "o", col = "#2b8cbe", lwd = 2)
	lines(quantiles[,3], x = newk, type = "o", col = "#74a9cf", lwd = 2)

	abline(v = newk, col = "lightgrey", lwd =1)
	abline(h = nullDistr[3], col = "#4daf4a", lwd = 3)
	abline(h = t0, col = "#e41a1c", lwd = 3)
	
	abline(v = x$uncertainty, col = "black")

	if(addMinimalStability){
		abline(a = nullDistr[2], b = a, lwd = 1, lty = 2)
		abline(a = nullDistr[1], b = a, lwd = 1, lty = 2)
	}

	# legend
	if(addLegend){
		tmp_t <- substitute(expression(t[tmp], "t*"), env = list(tmp = quant[3]))
		legend("bottomright",
			legend=c(paste(quant*100, "%-quantile", sep =""), eval(tmp_t)),
			lwd = 2,
			pch = c(1,1,1,NA,NA),
			col = c("#045a8d", "#2b8cbe", "#74a9cf", "#4daf4a", "#e41a1c"),
			bg = "white")
	}
}


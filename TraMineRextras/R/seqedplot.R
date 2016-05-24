seqedplot <- function(seqe, group=NULL, breaks=20, ages=NULL,title=NULL, type="survival", ignore=NULL,
	withlegend="auto",cex.legend=1, use.layout=(!is.null(group) | withlegend!=FALSE),legend.prop=NA, rows=NA, cols=NA, axes="all",
	xlab="time", ylab=ifelse(type=="survival", "survival probability", "mean number of events"), cpal=NULL, ...){
	
	if(!is.seqelist(seqe)) stop("seqe should be a seqelist. See help on seqecreate.")
	
	if(type=="survival" && all(seqelength(seqe)==-1)){
		stop(" [!] You should set observation time to use type='survival' (see seqelength)")
	}
	#if(is.null(ignore)) ignore <- character()
	## ==============================
	## Preparing if group is not null
	## ==============================
	if (!is.null(group)) {
		## Eliminate the unused levels
		group <- factor(group)
		nplot <- length(levels(group))
		gindex <- vector("list",nplot)
				
		for (s in 1:nplot)
			gindex[[s]] <- which(group==levels(group)[s])

		## Title of each plot
		if (!is.null(title))
			title <- paste(title,"-",levels(group))
		else
			title <- levels(group)
	}
	else {
		nplot <- 1
		gindex <- vector("list",1)
		gindex[[1]] <- 1:length(seqe)
	}
	
	## ===================
	## Defining the layout
	## ===================

	if (use.layout || !is.null(group) ) {
		## Saving graphical parameters
		savepar <- par(no.readonly = TRUE)
		on.exit(par(savepar))
		##lout <- TraMineR:::TraMineR.setlayout(nplot, rows, cols, withlegend, axes, legend.prop)
		lout <- TraMineRInternalLayout(nplot, rows, cols, withlegend, axes, legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)

		## Axis should be plotted or not ?
		xaxis <- 1:nplot==lout$axisp

		legpos <- lout$legpos
	}
	else {
		if(axes!=FALSE) xaxis <- TRUE
		else xaxis <- FALSE
		legpos <- NULL
	}
	agesmatrices <- list()
	allevents <- which(!levels(seqe) %in% ignore)
	numevent <- length(levels(seqe))
	nevent <- length(allevents)
	minage <- NA
	maxage <- NA
	for(event in 1:nevent){
        ##agematrix <- .Call("tmreventinseq", seqe, as.integer(allevents[event]), PACKAGE="TraMineR")
        agematrix <- TraMineRInternalSeqeage(seqe, as.integer(allevents[event]))
		agematrix[agematrix==-1] <- NA
		if (!is.null(ages)) {
			agematrix[agematrix<ages[1]] <- NA
			agematrix[agematrix>ages[2]] <- NA
		}
		suppressWarnings(minage <- min(c(minage, agematrix), na.rm=TRUE))
		suppressWarnings(maxage <- max(c(maxage, agematrix), na.rm=TRUE))
		#print(agematrix)
		agesmatrices[[event]] <- agematrix
	}
	if(is.null(ages))ages <- c(minage, maxage)
	if(is.null(cpal)) {
		if (nevent==2)
			cpal <- brewer.pal(3,"Dark2")[1:2]
		else if (nevent<=8)
			cpal <- brewer.pal(nevent,"Dark2")
		else if (nevent>8 & nevent<=12)
			cpal <- brewer.pal(nevent,"Set3")
		else
			cpal <- 1:nevent
	}
	## Finding correct cutpoints
	labs <- levels(cut(ages, breaks))
	cutpoints <- cbind(lower = as.numeric( sub("\\((.+),.*", "\\1", labs) ),
		upper = as.numeric( sub("[^,]*,([^]]*)\\]", "\\1", labs) ))
	cutpoints <- c(cutpoints[,1], max(cutpoints[,2]))
	onx <- cutpoints[-1]
	if (type=="hazard") {
		ony <- list()
		miny <- 1
		maxy <- 0
		for (np in 1:nplot) {
		ony[[np]] <- list()
			for (event in 1:nevent) {
				ony[[np]][[event]] <- as.numeric(xtabs(rep(seqeweight(seqe)[gindex[[np]]], ncol(agesmatrices[[event]]))~cut(agesmatrices[[event]][gindex[[np]],], breaks=cutpoints, ordered_result=TRUE))/length(gindex[[np]]))
				miny <- min(ony[[np]][[event]], miny)
				maxy <- max(ony[[np]][[event]], maxy)
			}
		}
	}
	for (np in 1:nplot) {
		if (type=="hazard") {
			doplot <- TRUE
			for (event in 1:nevent) {
				if (doplot) {
					plot(onx, ony[[np]][[event]], ylim=c(miny,maxy), col=cpal[event], type="l", xlab=xlab, ylab=ylab, main=title[np], ...)
					doplot <- FALSE
				}else{
					lines(onx, ony[[np]][[event]], col=cpal[event], ...)
				}
			}
		}
		else if (type=="survival") {
			## library(survival)
			doplot <- TRUE
			for (event in 1:nevent) {
				time <- agesmatrices[[event]][gindex[[np]],1]
				endobs <- seqelength(seqe)[gindex[[np]]]
				status <- !is.na(time)
				time[!status] <- endobs[!status]
				ony <- survfit(Surv(time, as.integer(status))~ 1)
				#ony <- ecdf(agesmatrices[[event]][gindex[[np]],1])
				if (doplot) {
					plot(ony, col=cpal[event], firstx=ages[1], xmax=ages[2], main=title[np], conf.int=FALSE, xlab=xlab, ylab=ylab, ...)
					doplot <- FALSE
				}else{
					lines(ony, col=cpal[event], ...)
				}
			}
		}
	}
	## Plotting the legend
	if (!is.null(legpos)) {
		## Extracting event names

		ltext <- levels(seqe)[allevents]

		##TraMineR:::TraMineR.legend(legpos, ltext, cpal, cex=cex.legend)
		TraMineRInternalLegend(legpos, ltext, cpal, cex=cex.legend)
	}

}

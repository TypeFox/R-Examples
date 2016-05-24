wcCmpCluster <- function(diss, weights=NULL, maxcluster, method="all", pam.combine=TRUE){
	if(maxcluster<2){
		stop(" [!] maxcluster should be greater than 2")
	}
	method <- tolower(method)
	hclustmethods <- c( "ward", "single", "complete", "average", "mcquitty", "median", "centroid")
	all.methods <- c(hclustmethods, "pam", "diana", "beta.flexible")
	noweights.methods <- c("diana", "beta.flexible")
	if(any(method=="all")){
		if(is.null(weights)){
			method <- all.methods
		} else{
			method <- c(hclustmethods, "pam")
		}
	}
	if(!is.null(weights) && any(method %in% noweights.methods)){
		stop(" [!] methods ", paste(noweights.methods, collapse=", "), " cannot be used with weights.")
	}
	ret <- list()
	if(any(method %in% hclustmethods)){
		dd <- as.dist(diss)
	}
	pamrange <- 2:maxcluster
	for(meth in method){
		if(meth !="pam"){
			if(meth %in% hclustmethods){
				hc <- hclust(dd, method=meth, members=weights)
			} else if(meth=="diana"){
				hc <- diana(diss, diss=TRUE)
			} else if(meth=="beta.flexible"){
				hc <- agnes(diss, diss=TRUE, method="flexible", par.method=0.625)
			}
			ret[[meth]] <- as.clustrange(hc, diss=diss, weights=weights, ncluster=maxcluster)
			if(pam.combine){
				ret[[paste("pam", meth, sep=".")]] <- wcKMedRange(diss, kvals=pamrange, weights=weights, initialclust=hc)
			}
		}else{
			ret[[meth]] <- wcKMedRange(diss, kvals=pamrange, weights=weights)
		}
	}
	allstats <- list()
	for(meth in names(ret)){
		allstats[[meth]] <- cbind(ret[[meth]]$stats, method=meth, ngroup=pamrange)
	}
	ret$param <- list(method=method, pam.combine=pam.combine, all.methods= names(ret), kvals=pamrange)
	ret$allstats <- do.call(rbind, allstats)
	class(ret) <- "clustrangefamily"
	return(ret)
}


print.clustrangefamily <- function(x, max.rank=1, ...){
	print(summary(x, max.rank=max.rank), ...)
}

summary.clustrangefamily <- function(object, max.rank=1, ...){
	statnames <- colnames(object$allstats)[1:(ncol(object$allstats)-2)]
	nstat <- length(statnames)
	clusterrank <- matrix(as.character(NA), nrow=nstat, ncol=max.rank*2)
	rownames(clusterrank) <- statnames
	clindices <- ((1:max.rank) *2) -1
	valindices <- ((1:max.rank) *2)
	cnames <- character(max.rank*2)
	cnames[clindices] <- paste(1:max.rank, "Cluster", sep=". ")
	cnames[valindices] <- paste(1:max.rank, " stat", sep=". ")
	colnames(clusterrank) <- cnames
	for(s in statnames){
		od <- order(object$allstats[, s], decreasing=(s!="HC"))[1:max.rank]
		clusterrank[s, clindices] <- paste(object$allstats[od, "method"], object$allstats[od, "ngroup"])
		clusterrank[s, valindices] <- object$allstats[od, s]
	}
	
	return(as.data.frame(clusterrank, check.names=FALSE))
}


crflegend <- function(pos, text, colors, cex = 1, ...){
	nbstat <- length(text)
    if (pos == "bottom") {
        if (nbstat > 6) 
            nbcol <- 3
        else nbcol <- 2
        leg.ncol <- ceiling(nbstat/nbcol)
    }
    else leg.ncol <- 1
    savepar <- par(mar = c(1, 1, 0.5, 1) + 0.1, xpd = FALSE)
    on.exit(par(savepar))
    plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
    legend(pos, legend = text, fill = colors, ncol = leg.ncol, 
        bty = "o", cex = cex, ...)

}

plot.clustrangefamily <- function(x, group="stat", method="all", pam.combine=FALSE, 
						stat="noCH", norm="none", 
						withlegend=TRUE, lwd=1, col=NULL, legend.prop=NA, rows=NA, 
						cols=NA, main=NULL, xlab="", ylab="", ...){
	savepar <- par(no.readonly = TRUE)
	on.exit(par(savepar))
	allstats <- colnames(x[[1]]$stats)
	if(length(stat)==1){
		if(stat=="all"){
			stat <- allstats
		}else if(stat=="noCH"){
			stat <- allstats[-grep("CH", allstats)]
		}
		else if(!stat %in%  c(allstats, "RHC")){
			stop(" [!] unknow statistic ", stat)
		}
	}else{
		stat <- stat[stat %in% c(allstats, "RHC")]
	}
	if(length(method)==1 && method=="all"){
		if(pam.combine){
		method <- x$param$all.methods
		} else{
			method <- x$param$method
		}
	} else if(length(method)==1 && method=="pam.combine"){
		if(!x$param$pam.combine){
			stop(" [!] you should have used pam.combine when computing the clustering.")
		}
		method <-  grep("pam", x$param$all.methods, value=TRUE)
	}
	intitle <- method
	if(group=="stat"){
		intitle <- stat
	}
	
	if(is.null(main)){
		main <- intitle
	} else if (length(main)==1 && length(intitle)>1){
		main <- paste(main, intitle, sep=" - ")
	} else if (length(main)!=length(intitle)){
		stop(" [!] You should provide one title (main argument) per plot.")
	}
	names(main) <- intitle
	
	if(group=="stat"){
		lout <- TraMineRInternalLayout(length(stat), rows, cols, withlegend, axes="all", legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)
		if(is.null(col)){
			if(length(method)>8){
				stop("[!] Too many clustering methods, you need to manually specify the colors (col argument).")
			}
			cols <- brewer.pal(length(method), "Dark2")
		} else {
			if(length(col) != length(method)){
				stop(" [!] You should specify at least one color per quality measure to plot.")
			}
			cols <- col
		}
		names(cols) <- method
		plotmat <- matrix(0, ncol=length(method), nrow=length(x$param$kvals))
		colnames(plotmat) <- method
		xlim <- range(x$param$kvals)
		for(st in stat){ ## each plots
			for(meth in method){
				if(st=="RHC"){
					plotmat[, meth] <- 1- x[[meth]]$stats[, "HC"]
				} else {
					plotmat[, meth] <- x[[meth]]$stats[, st]
				}
			}
			plotmat <- normalize.values.all(plotmat, norm)
			ylim <- range(unlist(plotmat), finite=TRUE)
			plot(xlim, ylim, xlim=xlim, ylim=ylim, type="n", xlab=xlab, ylab=ylab, main=main[st])
			for(meth in method){
				lines(x$param$kvals, plotmat[, meth], col=cols[meth], lwd=lwd, ...)
				##lines(x$param$kvals, plotmat[, meth], col=cols[meth], lwd=lwd)
			}
			
		}
		if(withlegend) {
			crflegend(lout$legpos, colors=cols, text=method)
		}
		
	} else { ## One plot per clustering method
		
		lout <- TraMineRInternalLayout(length(method), rows, cols, withlegend, axes="all", legend.prop)
	  	layout(lout$laymat, heights=lout$heights, widths=lout$widths)
		for(meth in method){
			plot(x[[meth]], withlegend=FALSE, lwd=lwd, col=col, stat=stat, norm=norm, main=main[meth], xlab=xlab, ylab=ylab, ...)
		}
		## Plot the legend
		
		if(is.null(col)){
			allnames <- allstats
			cols <- brewer.pal(length(allnames)+1, "Set3")[-2]
			names(cols) <- allnames
			cols["RHC"] <- cols["HC"]
			cols <- cols[stat]
		} else {
			if(length(col) != length(stat)){
				stop(" [!] You should specify at least one color per quality measure to plot.")
			}
			cols <- col
		}
		if(withlegend){
			crflegend(lout$legpos, colors=cols, text=stat)
		}
		
	}
}

# data(mvad)
# #Aggregating state sequence

##Creating state sequence object
# mvad.seq <- seqdef(mvad[, 17:86])

## COmpute distance using Hamming distance
# diss <- seqdist(biofam.seq, method="HAM")

##Ward clustering
# allClust <- clusterAllMethods(diss, weights=aggMvad$aggWeights, maxcluster=15)
# allClust <- clusterAllMethods(diss, maxcluster=15, method=c( "ward",  "average", "pam", "diana", "beta.flexible"))

# allClust
# ##Plot all statistics (standardized)
# plot(allClust, stat="all", norm="zscoremed", lwd=3)

# ##Plot HC, RHC and ASW
# plot(allClust, stat=c("HC", "RHC", "ASWw"), norm="zscore", lwd=3)
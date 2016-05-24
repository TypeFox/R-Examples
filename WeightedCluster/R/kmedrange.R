wcKMedRange <- function(diss, kvals, weights=NULL, R=1, samplesize=NULL, ...){
	if (inherits(diss, "dist")) {
		n <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		n <- nrow(diss)
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	ret <- list(R=R)
	ret$kvals <- kvals
	ret$clustering <- matrix(-1, nrow=n, ncol=length(kvals))
	ret$stats <-  matrix(-1, nrow=length(kvals), ncol=10)
	i <- 1
	for(k in kvals){
		cl <- wcKMedoids(diss=diss, k=k, cluster.only=TRUE, weights=weights, ...)
		ret$clustering[,i] <- cl
		i <- i+1
	}
	if(R==1){
		kendall <- .Call(wc_RClusterQualKendallFactory)
		for(i in 1:length(kvals)){
			stat <- wcClusterQualityInternal(diss=diss, clustering=ret$clustering[,i], weights=weights, kendall=kendall)
			ret$stats[i,] <- stat$stats
		}
		ret$stats <- as.data.frame(ret$stats)
		colnames(ret$stats) <- names(stat$stats)
	}else{
		ret$boot <- clustrangeboot(diss=diss, clustering=ret$clustering, weights=weights, R=R, samplesize=samplesize, simple=FALSE)
		for(i in 1:length(kvals)){
			ret$stats[i,] <- ret$boot[[i]][1, ]
		}
		ret$stats <- as.data.frame(ret$stats)
		colnames(ret$stats) <- colnames(ret$boot[[i]])
	}
	ret$clustering <- as.data.frame(ret$clustering)
	colnames(ret$clustering) <- paste("cluster", ret$kvals, sep="")
	rownames(ret$stats) <- paste("cluster", ret$kvals, sep="")
	class(ret) <- c("clustrange", class(ret))
	return (ret)
}

as.clustrange <- function(object, diss, weights=NULL, R=1, samplesize=NULL, ...){
	UseMethod("as.clustrange")
}

as.clustrange.hclust <- function(object, diss, weights=NULL, R=1, samplesize=NULL, ncluster=20, ...){
	
	if(ncluster<3){
		stop(" [!] ncluster should be greater than 2.")
	}
	if (is.null(n1 <- nrow(object$merge)) || n1 < 1) {
		stop("invalid 'object' (merge component)")
	}
    n <- n1 + 1
	if(ncluster > n){
		stop(" [!] ncluster should be less than ", n)
	}
	
	pred <- data.frame(Split2=factor(cutree(object, 2)))
	for(p in 3:ncluster){
		pred[, paste("Split", p, sep="")] <- factor(cutree(object, p))
	}
	object <- pred
	as.clustrange(object, diss=diss, weights=weights, R=R, samplesize=samplesize, ...)
}

as.clustrange.twins <- function(object, diss, weights=NULL, R=1, samplesize=NULL, ncluster=20, ...) {
	return(as.clustrange.hclust(object, diss=diss, weights=weights, ncluster=ncluster, R=R, samplesize=samplesize,...))
}

as.clustrange.default <- function(object, diss, weights=NULL, R=1,  samplesize=NULL,...){
	ret <- list()
	ret$clustering <- as.data.frame(object)
	numclust <- ncol(ret$clustering)
	ret$kvals <- numeric(numclust)
	ret$stats <-  matrix(-1, nrow=numclust, ncol=10)
	## print("BuildingKendall")
	kendall <- .Call(wc_RClusterQualKendallFactory)
	## print(class(kendall))
	## print(kendall)
	## print("Kendall")
	
	
	if(R==1){
		for(i in 1:numclust){
			ret$kvals[i] <- length(unique(ret$clustering[,i]))
			## print("Starting")
			cl <- wcClusterQualityInternal(diss, ret$clustering[,i], weights=weights, kendall=kendall)
			ret$stats[i,] <- cl$stats
		}
		ret$stats <- as.data.frame(ret$stats)
		colnames(ret$stats) <- names(cl$stats)
	}else{
		ret$boot <- clustrangeboot(diss=diss, clustering=ret$clustering, weights=weights, R=R, samplesize=samplesize, simple=FALSE)
		ret$meant <- ret$stats
		ret$stderr <- ret$stats
		for(i in 1:numclust){
			ret$kvals[i] <- length(unique(ret$clustering[,i]))
			ret$stats[i,] <- ret$boot[[i]][1, ]
			ret$meant[i,] <- colMeans(ret$boot[[i]])
			ret$stderr[i,] <- apply(ret$boot[[i]], 2L, function(x) sqrt(var(x)))
		}
		ret$stats <- as.data.frame(ret$stats)
		ret$meant <- as.data.frame(ret$meant)
		ret$stderr <- as.data.frame(ret$stderr)
		colnames(ret$stats) <- colnames(ret$boot[[i]])
		colnames(ret$meant) <- colnames(ret$boot[[i]])
		colnames(ret$stderr) <- colnames(ret$boot[[i]])
		rownames(ret$meant) <- paste("cluster", ret$kvals, sep="")
		rownames(ret$stderr) <- paste("cluster", ret$kvals, sep="")
	}
	colnames(ret$clustering) <- paste("cluster", ret$kvals, sep="")
	rownames(ret$stats) <- paste("cluster", ret$kvals, sep="")
	class(ret) <- c("clustrange", class(ret))
	return(ret)
}
print.clustrange <- function(x, digits=2, bootstat=c("t0", "mean", "stderr"), ...){
	if(!is.null(x$boot) && bootstat!="t0"){
		if(bootstat=="mean"){
			x <- round(x$meant, digits)
		}else{
			x <- round(x$stderr, digits)
		}
	} else {
		x <- round(x$stats, digits)
	}
	print(x, ...)
	
}

normalize.values.all <- function(stats, norm){
	for(i in 1:ncol(stats)){
		stats[,i] <- normalize.values(stats[, i], norm)
	}
	return(stats)
}
normalize.values <- function(stats, norm){
	if(norm == "range") {
		stats <- (stats-min(stats))/(max(stats)-min(stats))
	} else if (norm=="zscore") {
		stats <- (stats - mean(stats ))/(sqrt(var(stats)))
	}else if (norm=="zscoremed") {
		stats <- (stats - median(stats))/(median(abs(stats - median(stats))))
	}
	return(stats)
}
normalize.values.matrix <- function(stats, norm){
	st <- normalize.values(c(stats), norm)
	dim(st) <- dim(stats)
	return(st)
} 

plot.clustrange <- function(x, stat="noCH", legendpos="bottomright", norm="none", 
							withlegend=TRUE, lwd=1, col=NULL, ylab="Indicators", 
							xlab="N clusters", conf.int=0.9, ci.method="none", ci.alpha=.3, line="t0", ...){
	kvals <- x$kvals
	if(length(stat)==1){
		if(stat=="all"){
			stats <- x$stats
			if(norm=="none"){
				norm <- "range"
			}
		}else if(stat=="noCH"){
			stats <- x$stats[,-grep("CH", colnames(x$stats))]
		}
		else{
			stats <- x$stats[,stat]
		}
	}else{
		if("RHC" %in% stat){
			stats <- x$stats[ , stat[stat != "RHC"]]
			stats[, "RHC"] <- 1 - x$stats[ , "HC"]
		}
		else{
			stats <- x$stats[, stat]
		}
	}
	
	if(is.null(col)){
		allnames <- colnames(x$stats)
		cols <- brewer.pal(length(allnames)+1, "Set3")[-2]
		names(cols) <- allnames
		cols["RHC"] <- cols["HC"]
		cols <- cols[colnames(stats)]
	} else {
		if(length(col) < ncol(stats)){
			stop(" [!] You should specify at least one color per quality measure to plot.")
		}
		cols <- col
	}
	plot.ci <-  ci.method!="none" && !is.null(x$boot)
	if(plot.ci) {
		upper <- stats
		lower <- stats
		confvec <- c((1-conf.int)/2, 1-(1-conf.int)/2)
		getline <- function(x){
			if(line=="t0"){
				return(x[1])
			}else if(line=="mean"){
				return(mean(x))
			}else if(line=="median"){
				return(median(x))
			}
		}
		borne <- list()
		for(l in colnames(stats)){
			st <- normalize.values.matrix(sapply(x$boot, function(x)x[ , l]), norm=norm)
			stats[, l] <- apply(st, 2, getline)
			if(ci.method=="norm"){
				borne[[l]] <- apply(st, 2, function(x) {mean(x)+qnorm(confvec)*sqrt(var(x))})
			}
			else if(ci.method=="perc"){
				borne[[l]] <- apply(st, 2, function(x) quantile(x, confvec))
			}
		}
		ylim <- range(unlist(c(stats, borne)), finite=TRUE)
	}else{
		stats <- normalize.values.all(stats, norm)
		ylim <- range(unlist(stats), finite=TRUE)
	}
	plot(kvals, xlim=range(kvals, finite=TRUE), ylim=ylim, type="n", ylab=ylab, xlab=xlab, ...)
	labels <- paste(colnames(stats), "(", round(apply(stats, 2, min), 2),"/", round(apply(stats, 2, max),2), ")")	
	names(labels) <- colnames(stats)
	for(l in colnames(stats)){
		ss <- stats[,l]
		lines(kvals, ss, col=cols[l], lwd=lwd, ...)
		if(plot.ci){
			polygon(c(rev(kvals), kvals), c(rev(borne[[l]][1, ] ), borne[[l]][2, ] ), col = adjustcolor(cols[l], alpha.f=ci.alpha), border = NA)
		}
		
	}
	if(withlegend) {
		legend(legendpos, fill=cols[1:ncol(stats)], legend=labels)
	}
}

summary.clustrange <- function(object, max.rank=1, ...){
	nstat <- ncol(object$stats)
	clusterrank <- matrix(NA, nrow=nstat, ncol=max.rank*2)
	rownames(clusterrank) <- colnames(object$stats)
	clindices <- ((1:max.rank) *2) -1
	valindices <- ((1:max.rank) *2)
	cnames <- character(max.rank*2)
	cnames[clindices] <- paste(1:max.rank, "N groups", sep=". ")
	cnames[valindices] <- paste(1:max.rank, " stat", sep=". ")
	colnames(clusterrank) <- cnames
	for(s in colnames(object$stats)){
		od <- order(object$stats[, s], decreasing=(s!="HC"))[1:max.rank]
		clusterrank[s, clindices] <- object$kvals[od]
		clusterrank[s, valindices] <- object$stats[od, s]
	}
	
	return(as.data.frame(clusterrank, check.names=FALSE))
}

wcClusterQuality <- function(diss, clustering, weights=NULL){
	return(wcClusterQualityInternal(diss=diss, clustering=clustering, weights=weights, kendall=NULL))
}
wcClusterQualityInternal <- function(diss, clustering, weights=NULL, kendall=NULL){
	if (inherits(diss, "dist")) {
		isdist <- TRUE
		nelements <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		isdist <- FALSE
		nelements <- nrow(diss)
		if(ncol(diss)!=nelements){
			stop("[!] diss should be a squared matrix or a dist object.")
		}
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	clusterF <- factor(clustering)
	if(nlevels(clusterF) < 2){
		stop("[!] The clustering should have at least two different values.")
	}
	clustering <- as.integer(as.integer(clusterF)-1)
	if(length(clustering)!=nelements|| length(weights)!=nelements){
		stop("[!] different number of elements in diss, clustering and/or weights arguments.")
	}
	
	ncluster <- max(clustering)+1
	if(is.null(kendall)){
		cq <- .Call(wc_RClusterQual, diss, clustering, as.double(weights), as.integer(ncluster), as.integer(isdist), as.integer(0))
	}else{
		cq <- .Call(wc_RClusterQualKendall, diss, clustering, as.double(weights), as.integer(ncluster), as.integer(isdist), kendall)
	}
	names(cq[[1]]) <-c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", "R2sq", "HC")
	dim(cq[[2]]) <- c(nlevels(clusterF), 2)
	
	rownames(cq[[2]]) <-levels(clusterF)
	colnames(cq[[2]]) <- c("ASW", "ASWw")
	names(cq) <- c("stats", "ASW")
	if(any(xtabs(weights~clustering)<1)){
		cq$ASW[, "ASW"] <- NA
		cq$stats["ASW"] <- NA
		warning(" [!] ASW can not be computed because at least one cluster has less than one observation.\n")
	}
	return(cq)
#define ClusterQualHPG 0
#define ClusterQualHG 1
#define ClusterQualHGSD 2
#define ClusterQualASWi 3
#define ClusterQualASWw 4
#define ClusterQualF 5
#define ClusterQualR 6
#define ClusterQualF2 7
#define ClusterQualR2 8
#define ClusterQualHC 9
#define ClusterQualNumStat 10

}

wcSilhouetteObs <- function(diss, clustering, weights=NULL, measure="ASW"){
	if(!all(measure %in% c("ASW", "ASWw"))){
		stop(" [!] Unknow silhouette measure. It should be one of 'ASW', 'ASWw'.")
	}
	if (inherits(diss, "dist")) {
		isdist <- TRUE
		nelements <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		isdist <- FALSE
		nelements <- nrow(diss)
		if(ncol(diss)!=nelements){
			stop("[!] diss should be a squared matrix or a dist object.")
		}
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	clusterF <- factor(clustering)
	if(nlevels(clusterF) < 2){
		stop("[!] The clustering should have at least two different values.")
	}
	clustering <- as.integer(as.integer(clusterF)-1)
	if(length(clustering)!=nelements|| length(weights)!=nelements){
		stop("[!] different number of elements in diss, clustering and/or weights arguments.")
	}
	if(any(measure=="ASW") && any(xtabs(weights~clustering)<1)){
		if(length(measure)==1){
			stop(" [!] ASW can not be computed because at least one cluster has less than one observation. Use measure='ASWw'.")
		}
		warning(" [!] ASW can not be computed because at least one cluster has less than one observation. Only 'ASWw' is computed.")
		measure <- "ASWw"
	}
	ncluster <- max(clustering)+1
	ret <- .Call(wc_RClusterComputeIndivASW, diss, clustering, as.double(weights), as.integer(ncluster), as.integer(isdist))
	asw <- data.frame(ASW=ret[[1]], ASWw=ret[[2]])
	return(asw[, measure])
	#define ClusterQualHPG 0
#define ClusterQualHG 1
#define ClusterQualHGSD 2
#define ClusterQualASW 3
#define ClusterQualF 4
#define ClusterQualR 5
#define ClusterQualF2 6
#define ClusterQualR2 7
}



clustrangeboot<- function(diss, clustering, weights=NULL, R=999, samplesize=NULL, simple=FALSE){
	if(simple){
		statname <- c("PBC", "CH", "R2", "CHsq", "R2sq")
	}else{
		statname <- c("PBC", "HG", "HGSD", "ASW", "ASWw", "CH", "R2", "CHsq", "R2sq", "HC")
	}
	nstat <- length(statname)
	if (inherits(diss, "dist")) {
		isdist <- TRUE
		nelements <- attr(diss, "Size")
	}else if(is.matrix(diss)){
		isdist <- FALSE
		nelements <- nrow(diss)
		if(ncol(diss)!=nelements){
			stop("[!] diss should be a squared matrix or a dist object.")
		}
	} else {
		stop("[!] diss should be a squared matrix or a dist object.")
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	if(!is.data.frame(clustering)){
		clustering <- data.frame(clustering)
	}
	if(nrow(clustering)!=nelements|| length(weights)!=nelements){
		stop("[!] different number of elements in diss, clustering and/or weights arguments.")
	}
	clustmat <- matrix(0L, nrow=nrow(clustering), ncol=ncol(clustering))
	nclusters <- integer(ncol(clustering))
	ans <- vector("list", length=ncol(clustering))
	for(i in 1:ncol(clustering)){
		ans[[i]] <- matrix(0.0, ncol=nstat, nrow=R+1, dimnames=list(NULL, statname))
		clusterF <- factor(clustering[, i])
		if(nlevels(clusterF) < 2){
			stop("[!] Each clustering should have at least two different values.")
		}
		clustmat[, i] <- as.integer(as.integer(clusterF)-1)
		nclusters[i] <- as.integer(nlevels(clusterF))
	}

	
	totweights <- sum(weights)
	prob <- weights/totweights
	if(is.null(samplesize)){
		samplesize <- as.integer(floor(totweights))
	}
	
	internalsample <- function(){
		return(as.integer(sample.int(nelements, size=samplesize, replace=TRUE, prob=prob)-1L))
	}
	
	bts <- .Call("RClusterQualBootSeveral", ans, diss, clustmat, as.double(weights), 
										as.integer(nclusters), as.integer(R+1),  quote(internalsample()), 
										environment(), as.integer(samplesize), as.integer(isdist), as.integer(simple), PACKAGE="WeightedCluster")
	
	class(ans) <- "clustrangeboot"
	return(ans)
}

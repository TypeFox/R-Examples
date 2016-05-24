wcKMedoids <- function(diss, k, weights=NULL, npass = 1, initialclust=NULL, method="PAMonce", cluster.only = FALSE, debuglevel=0) {
	if(is.character(method) && !is.na(pmatch(tolower(method), c("kmedoids", "pam", "pamonce")))){
		method <- pmatch(tolower(method), c("kmedoids", "pam", "pamonce"))
	}
	if(!(is.numeric(method) && method %in% 1:3)){
		stop(" [!] unknow clustering method ", method)
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
	InternalRandomSample <- function(){
		return(as.integer(sample.int(nelements, k)-1))
	}
	if(is.null(weights)){
		weights <- as.double(rep(1.0, nelements))
	}
	if(length(weights)!=nelements){
		stop("[!] weights should be a vector of length ", nelements)
	}
	if(is.null(initialclust)) {
		initialclust <- InternalRandomSample()
	} else {
		if(inherits(initialclust, "hclust") || inherits(initialclust, "twins")){
			initialclust <- cutree(initialclust, k)
		}
		if(length(initialclust)==nelements){
			initialclust <- disscenter(diss, group=initialclust, weights=weights, medoids.index="first")
			if(length(initialclust) != k ){
				stop(" [!] initialclust should be a vector of cluster membership with k=", k, " different groups.")
			}
		}
		initialclust <- initialclust -1
		npass <- 0
	}
	if(length(initialclust) != k ){
		stop(" [!] initialclust should be a vector of medoids index of length :", k)
	}
	if(any(initialclust>=nelements | initialclust<0)){
		stop(" [!] starting medoids should be in 1:", nelements)
	}
	if(npass < 0) {
		stop(" [!] npass should be greater than 0")
	}
	if(k < 2 || k > nelements) {
		stop(" [!] k should be in [2,",nelements,"]")
	}
	
	ret <- .Call(wc_RKmedoids, as.integer(nelements),  diss, quote(InternalRandomSample()), environment(), as.integer(initialclust), as.integer(npass), as.double(weights), as.integer(method), as.integer(debuglevel), as.integer(isdist))
	names(ret) <- c("clustering", "info")
	## Taking care of C style indices
	ret$clustering <- ret$clustering+1
	if(cluster.only){
		return(ret$clustering)
	}
	allstat <- wcClusterQuality(diss=diss, clustering=ret$clustering, weights=weights)
	ret$stats <- allstat$stats
	ret$ASW <- allstat$ASW
	names(ret$info) <- c("Sum of Distances", "N. Found", "Method")
	class(ret) <- c("kmedoids", class(ret))
	return(ret)
}

print.kmedoids <- function(x, ...){
	##cat("Basic info:\n")
	##print(as.data.frame(t(x$info)), row.names = FALSE, ...)
	cat("Cluster quality:\n")
	print(as.data.frame(t(x$stats)), row.names = FALSE, ...)
	cat("\nPer cluster Weighted Average Silhouette Width:\n")
	print(as.data.frame(t(x$ASW)), row.names = FALSE, ...)
}


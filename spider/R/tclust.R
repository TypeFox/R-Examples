tclust <- function(distobj, threshold = 0.01){
	distobj <- as.matrix(distobj)
	out <- unique(apply(distobj, 2, function(x) which(x < threshold)))
	for(i in 1:length(out)){
	for(j in 1:length(out)){
		if(length(intersect(out[[i]], out[[j]])) > 0) out[[i]] <- union(out[[i]], out[[j]]) else out[[i]] <- out[[i]]
		}
	}
	unique(lapply(out, sort))
}
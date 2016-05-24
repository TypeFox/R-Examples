get.cluster <-
function(A, iso.seq, criteria.list, clusters.list){
	
	n = dim(A)[1]
	
	ratio.list = criteria.list$ratio.list
	normalised.list = criteria.list$normalised.list
	
	ratio.location = which.min(ratio.list)
	ratio.cluster = clusters.list[, ratio.location]
	
	normalised.location = which.min(normalised.list)
	normalised.cluster = clusters.list[, normalised.location]
	
	
	if(length(iso.seq) != 0) noniso.seq = (1:n)[-iso.seq]
	else noniso.seq = (1:n)
	
	final.ratio.cluster = rep(0, dim(A)[1])
	final.ratio.cluster[noniso.seq] = ratio.cluster
	
	final.normalised.cluster = rep(0, dim(A)[1])
	final.normalised.cluster[noniso.seq] = normalised.cluster
	
	
	return(list(final.ratio.cluster = final.ratio.cluster, ratio.location = ratio.location, final.normalised.cluster = final.normalised.cluster, normalised.location = normalised.location))

}

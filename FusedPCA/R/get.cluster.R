get.cluster <- function(A, iso.seq, cut.list, clusters.list, mod.list){
	
	n = dim(A)[1]
	
	ratio.list = cut.list$ratio.list
	normalised.list = cut.list$normalised.list
	
	ratio.location = which.min(ratio.list)
	ratio.cluster = clusters.list[, ratio.location]
	
	normalised.location = which.min(normalised.list)
	normalised.cluster = clusters.list[, normalised.location]
	
	mod.dcbm.list = mod.list$mod.dcbm.list
	mod.sbm.list = mod.list$mod.sbm.list
	
	mod.dcbm.location = which.min(mod.dcbm.list)
	mod.dcbm.cluster = clusters.list[, mod.dcbm.location]

	mod.sbm.location = which.min(mod.sbm.list)
	mod.sbm.cluster = clusters.list[, mod.sbm.location]
	
	noniso.seq = (1:n)[-iso.seq]
	
	final.ratio.cluster = rep(0, dim(A)[1])
	final.ratio.cluster[noniso.seq] = ratio.cluster
	
	final.normalised.cluster = rep(0, dim(A)[1])
	final.normalised.cluster[noniso.seq] = normalised.cluster
	
	final.mod.dcbm.cluster = rep(0, dim(A)[1])
	final.mod.dcbm.cluster[noniso.seq] = mod.dcbm.cluster

	final.mod.sbm.cluster = rep(0, dim(A)[1])
	final.mod.sbm.cluster[noniso.seq] = mod.sbm.cluster
	
	return(list(final.ratio.cluster = final.ratio.cluster, ratio.location = ratio.location, final.normalised.cluster = final.normalised.cluster, normalised.location = normalised.location, final.mod.dcbm.cluster = final.mod.dcbm.cluster, mod.dcbm.location = mod.dcbm.location, final.mod.sbm.cluster = final.mod.sbm.cluster, mod.sbm.location = mod.sbm.location))

}
fpca <- function(A, maxsteps = 200, tol = 1e-3, normalised = T, K = 2, score = F, ridge = T, approx = F){
	
	temp.obj = fpca.start(A = A, maxsteps = maxsteps, tol = tol, normalised = normalised, K = K, score = score, ridge = ridge, approx = approx)
	temp.cluster = fpca.cluster(obj = temp.obj, K = K, score = score)
	temp.cut = fpca.cut(A = A, obj = temp.obj, fpca.cluster = temp.cluster, K = K)
	temp.mod = fpca.mod(A = A, obj = temp.obj, fpca.cluster = temp.cluster, K = K)
	temp = get.cluster(A = A, iso.seq = temp.obj$iso.seq, cut.list = temp.cut, clusters.list = temp.cluster, mod.list = temp.mod)
	
	return(list(final.ratio.cluster = temp$final.ratio.cluster, ratio.location = temp$ratio.location, final.normalised.cluster = temp$final.normalised.cluster, normalised.location = temp$normalised.location, final.mod.dcbm.cluster = temp$final.mod.dcbm.cluster, mod.dcbm.location = temp$mod.dcbm.location, final.mod.sbm.cluster = temp$final.mod.sbm.cluster, mod.sbm.location = temp$mod.sbm.location))
	
}
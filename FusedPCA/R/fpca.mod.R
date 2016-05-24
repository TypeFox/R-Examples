fpca.mod <- function(A, obj, fpca.cluster, K = 2, iso.seq){
	
	n.loc = dim(fpca.cluster)[2]
	mod.dcbm.list = c()
	mod.sbm.list = c()
	
	n = dim(A)[1]
	if(length(obj$iso.seq) > 0) noniso = (1:n)[-obj$iso.seq]
	else noniso = (1:n)
	A.noniso = A[noniso, noniso]
	for(i in 1: n.loc){
		temp = single.mod(A.noniso, clusters = fpca.cluster[, i], K = K)
		mod.dcbm.list[i] = temp$mod.dcbm
		mod.sbm.list[i] = temp$mod.sbm
	}
	return(list(mod.dcbm.list = mod.dcbm.list, mod.sbm.list = mod.sbm.list))
}
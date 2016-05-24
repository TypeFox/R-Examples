fpca.cut <- function(A, obj, fpca.cluster, K = 2, iso.seq){
	
	n.loc = dim(fpca.cluster)[2]
	ratio.list = c()
	normalised.list = c()

	n = dim(A)[1]
	if(length(obj$iso.seq) > 0) noniso = (1:n)[-obj$iso.seq]
	else noniso = (1:n)
	A.noniso = A[noniso, noniso]
	for(i in 1: n.loc){
		temp = single.cut(A.noniso, clusters = fpca.cluster[, i], K = K)
		ratio.list[i] = temp$ratio.count
		normalised.list[i] = temp$normalised.count
	}
	return(list(ratio.list = ratio.list, normalised.list = normalised.list))
}

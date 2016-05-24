fcd.criteria <-
function(A, fcd.cluster, K = 2, iso.seq){
	
	### preparation
	n.loc = dim(fcd.cluster)[2]
	ratio.list = c()
	normalised.list = c()
	
	### main
	n = dim(A)[1]
	if(length(iso.seq) > 0) noniso = (1:n)[-iso.seq]
	else noniso = (1:n)
	A.noniso = A[noniso, noniso]
	for(i in 1: n.loc){
		temp1 = single.cut(A.noniso, clusters = fcd.cluster[, i], K = K)
		ratio.list[i] = temp1$ratio.count
		normalised.list[i] = temp1$normalised.count
		
	}
	
	### return
	return(list(ratio.list = ratio.list, normalised.list = normalised.list))
}

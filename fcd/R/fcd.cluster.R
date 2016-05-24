fcd.cluster <-
function(obj, K = 2){
	
	### preparation
	beta.combind = obj
	min.length = length(beta.combind)
	n = dim(beta.combind[[1]])[1]
	cluster.list = matrix(0, nrow = n, ncol = min.length)
	
	### clustering
	for(i in 1: min.length){
		temp1 = unique(beta.combind[[i]], margin = 2)
		if(dim(temp1)[1] >= K){
			k.means = kmeans(beta.combind[[i]], centers = K)
			cluster.list[, i] = k.means$cluster
		}
	}

	### return
	return(cluster.list)
}

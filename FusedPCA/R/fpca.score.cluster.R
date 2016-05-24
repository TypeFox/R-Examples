fpca.score.cluster <- function(obj, K = 2){
	
	n = dim(obj)[1]
	n.loc = dim(obj)[3]
	cluster.list = matrix(0, nrow = n, ncol = n.loc)

	for(i in 1: n.loc){
		
		temp1 = unique(as.matrix(obj[, , i]), margin = 2)
		if(dim(temp1)[1] >= K){
			k.means = kmeans(obj[, , i], centers = K)
			cluster.list[, i] = k.means$cluster
		}
	}
	
	return(cluster.list)
	
}
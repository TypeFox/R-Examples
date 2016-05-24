single.mod <- function(A, clusters, K = 2){
	
	eps = 1e-5
	n.list = c()
	O.vector = c()
	O.matrix = matrix(0, nrow = K, ncol = K)
	for(i in 1:K){
		label.1 = which(clusters == i)
		n.list[i] = length(label.1)+eps
		for(j in 1:K){
			label.2 = which(clusters == j)
			O.matrix[i, j] = sum(A[label.1, label.2])+eps
		}
		O.vector[i] = sum(A[label.1, ])+eps
	}
	
	mod.dcbm = 0
	mod.sbm = 0
	for(i in 1:K){
		for(j in 1:K){
			mod.dcbm = mod.dcbm + O.matrix[i,j]*log(O.matrix[i,j]/(O.vector[i] * O.vector[j]))
			mod.sbm = mod.sbm + O.matrix[i,j]*log(O.matrix[i,j]/(n.list[i] * n.list[j]))
		}
	}
	
	return(list(mod.dcbm = mod.dcbm, mod.sbm = mod.sbm))	
}
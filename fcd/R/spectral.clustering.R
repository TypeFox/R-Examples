spectral.clustering <- function(A, normalised = TRUE, score = FALSE, K = 2, adj = FALSE){

	### preparation 
	n = dim(A)[1]
	iso.A = isolate(A)
	iso.seq = iso.A$isolate
	noniso.seq = iso.A$nonisolate
	A.noniso = A[noniso.seq, noniso.seq]
	labels = rep(0, n)
	
	### svd
	n = dim(A.noniso)[1]
	eig.A = eigen(A.noniso)
	if(score == F){
		U = matrix(0, nrow = n, ncol = K)	
	}
	if(score == T){
		U = matrix(0, nrow = n, ncol = (K-1))
	}
	
	### get U matrix
	if(adj == F & score == F){
		L = laplacian(A = A.noniso, normalised = normalised)
		eig = eigen(L, symmetric = TRUE)
	
		for(i in 1:K){
			U[, i] = eig$vectors[, (n + 1 - i)]
		}
	}
	
	if(adj == T & score == F){
		ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
		for(i in 1:K){
			U[, i] = ordered.vec[, (n + 1 - i)]
		}
	}
	
	if(score == T){
		ordered.vec = eig.A$vectors[, order(abs(eig.A$values), decreasing = T)]
		benchmark = ordered.vec[,1] + 1e-5
		for(i in 2:K){
			U[, (i-1)] = eig.A$vectors[, order(abs(eig.A$values), decreasing = F)][, i]/benchmark
		}
	}
	
	### k-means
	U = scale(U, center = F)
	temp = unique(U, margin = 2)
	if(dim(temp)[1] < K){stop('FAIL!')}
	k.means = kmeans(U, centers = K)
	labels[noniso.seq] = k.means$cluster
	
	return(labels)
}




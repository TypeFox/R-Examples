find_nn_eps <-
function(X,eps){
	
	#calculate distance matrix
	nns <- as.matrix(dist(X))

	#choose neighbours using eps environment
	nns <- nns<eps
	diag(nns) <- FALSE
	
	return(nns)
}


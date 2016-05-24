DirichNormSat <-
function(nvec,bvec){
	if (length(nvec) != length(bvec) ) stop("DirichNormSat: Dimension of nvec and bvec differ\n") 
	DirichNormSat <- lfactorial(sum(nvec)) - sum(lfactorial(nvec)) + 
		         lgamma(sum(bvec)) - sum(lgamma(bvec)) + 
			 sum(lgamma(nvec+bvec)) - lgamma(sum(nvec)+sum(bvec))
	DirichNormSat <- exp(DirichNormSat)
	DirichNormSat
}


## ========================================
## Read a TDA output into a distance matrix
## ========================================

read.tda.mdist <- function(file) {
	tmp <- read.table(file=file,header=FALSE)
	Y <- nrow(tmp)
	d <- 1+8*Y
	n <- (1+sqrt(d))/2

	cat(" [>] reading a",n,"x",n,"matrix\n")
	m <- matrix(nrow=n,ncol=n)
	diag(m) <- 0

	for (i in 1:Y) {
		m[tmp[i,"V1"],tmp[i,"V2"]] <- tmp[i,"V5"]
	}

	m[upper.tri(m)] <- t(m)[upper.tri(m)]

	return(m)
}
		
			

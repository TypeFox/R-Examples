# 
# compute the stationary distribution of a transition matrix
# 
# Ingmar Visser, 09-2009
# 

stationary <- function(tpm) {
	if(!(nrow(tpm)==ncol(tpm))) stop("stationary distribution only defined for square matrices.")
	nr <- nrow(tpm)
	rs <- all.equal(rowSums(tpm),rep(1,nr))
	if(!(rs==TRUE)) stop("rows of the transition matrix should sum to unity")
	e1 <- as.double(eigen(t(tpm))$vectors[,1])
	e1 <- e1/sum(e1)
	return(e1)
}

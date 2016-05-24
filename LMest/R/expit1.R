expit1 <- function(lp,ref=1){
	
# compute expit with respect to category ref
	k = length(lp)+1
	G = matrix(diag(k)[,-ref],k,k-1)
	p = exp(G%*%lp); p = p/sum(p)
	p = as.vector(p)
	Der = (diag(p)-p%o%p)%*%G
	out = list(p=p,Der=Der)
		
}
lginitc <-
function(p,q,G,sgc){
	d <- matrix(0, nrow=p, ncol=p)
	lg <- array(0, dim=c(p, q, G))
	eigval <- eigen(sgc,symmetric=TRUE)$values
			if(min(eigval)<0){
				#message("A NEGATIVE EIGENVALUE!? 
				#I've taken the absolute value instead..."); 
				eigval<-abs(eigval)
			}
	diag(d) <- sqrt(eigval)
	fullLg <- eigen(sgc,symmetric=TRUE)$vectors %*% d
	dum <- fullLg[,1:q]
	for(g in 1:G){
		lg[,,g] <- dum
	}
	lg
}

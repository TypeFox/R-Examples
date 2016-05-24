lginitu <-
function(p,q,G,sg){
	d <- matrix(0,p,p)
	lg <- array(0, dim=c(p, q, G))
	for(g in 1:G){	
		eigval <- eigen(sg[,,g],symmetric=TRUE)$values
			if(min(eigval)<0){
				#message("A NEGATIVE EIGENVALUE!? 
				#I've taken the absolute value instead..."); 
				eigval<-abs(eigval)
			}
		diag(d) <- sqrt(eigval)
		fullLg <- eigen(sg[,,g],symmetric=TRUE)$vectors %*% d
		lg[,,g] <- fullLg[,1:q]
	}
	lg
}

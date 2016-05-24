"GARI" <- 
function(ADM,EADM){
	n <- nrow(ADM)
	N <- n-1
	tmp <- rep(0,n)
	for(i in 1:n){
		k <- sum(ADM[i,])
		tmp[i] <- N + 2*k*(k-N)/N
	}
	ECN <- sum(tmp)
	X <- sum(ADM==EADM) - n
	gari <- (X- ECN )/( n*(n-1) - ECN )
	return( gari )
}
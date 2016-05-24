

deriva_grad <- function(b,grad){
	m <- length(b)
	bh2 <- b
	bh <- b 
	v <- rep(0,m*(m+1)/2)
	k <- 0	
	for(j in 1:m){
		for(i in 1:j){
			k <- k+1
			bh <- b 
			bh2 <- b
			thj <- -max(1e-7, 1e-4 * abs(b[j]))
			bh2[j] <- b[j]-thj
			bh[j] <- b[j]+thj
			v[k] <- (grad(bh)[i]-grad(bh2)[i])/(2*(-thj))
		}
	}	
	result <- list(hessian=v)
	return(result)
}


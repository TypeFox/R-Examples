sim.michael = function(N){
	stopifnot(N>2)
	P = 0.25
	Q = 0.25
	n <- min(N,510)
	
	for(n in 2:n){
		P[n] = (choose(2*n,n)/4^n)^2
 		Q[n] = P[n] - sum(Q[1:(n-1)]*P[(n-1):1])
 	}
 	b <- sum(Q[1:n] * log(c(1:n))) - pi*log(log(n))
 	if(N > n){
		bsp <- .Call("simmich",P,Q,as.integer(n) ,as.integer(N) )
		print(bsp[2])
		b <- bsp[1]
	}
	
	return(b)
}





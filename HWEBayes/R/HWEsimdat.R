HWEsimdat <-
function(npop,q,f){
	 pmin <- min(q)
         fmin <- -pmin/(1-pmin)  
	 if (f<fmin | f >1) stop("HWsimdat: Illegal value of f")
	 k <- length(q)
	 ncell <- k*(k+1)/2
	 p <- rep(0,ncell)
	 count <- 1
	 for (i in 1:k){
	     for (j in i:k){
	     	 if (i==j) {p[count] <- q[i]^2 + f*q[i]*(1-q[i]); 
		    	    count <- count+1}
	         if (j>i){p[count] <- 2*q[i]*q[j]*(1-f); 
		    	    count <- count+1}
	     }
	 }
	 nvec <- rmultinom(n=1, size=npop, prob=p)
	 list(nvec=nvec)
}


DirichSampHWE <-
function(nvec,bvec0,nsim){
# first convert to allele counts mvec
	   k <- .5*(-1+sqrt(1+8*length(nvec)))
	   if (k != length(bvec0) ) stop("DirichSampHWE: Dimension of nvec and bvec differ\n") 
# turn nvec into kxk matrix
  	   nmat <- matrix(0,nrow=k,ncol=k); count <- 1
	   mvec <- rep(0,k)
	   for(i in 1:k){
      	   	 for (j in 1:k){
      	  	     if (i<=j){
      	     	     	nmat[i,j] <- nvec[count]
	     	     	count <- count+1}}}
		 for(i in 1:k){
      	  	   for (j in 1:k){
      	  	     if (i>j) nmat[i,j] <- nmat[j,i]
		   }
		 }
		 for (i in 1:k){
		     mvec[i] <- 2*nmat[i,i]
		     for (j in 1:k){
		     	 if (j != i) mvec[i] <- mvec[i] + nmat[i,j]}}
#cat("Counts: ",nmat,"\n")
	#   cat("Allele counts: ",mvec,"\n")
# posterior samples for the allele probs
    	   pvec <-  rdirichlet(nsim,mvec+bvec0)
	   list(pvec=pvec)
}


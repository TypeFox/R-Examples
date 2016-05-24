DirichNormHWE <-
function(nvec,bvec0){
	 k <- .5*(-1+sqrt(1+8*length(nvec)))
	 if (length(bvec0) != k) stop("DirichNormHWE: Dimension of bvec0 not equal to k\n")
	 nsum <- sum(nvec)
	 bsum <- sum(bvec0)
	 DirichNormHWE <- lfactorial(nsum) - sum(lfactorial(nvec)) + lgamma(bsum)  - sum(lgamma(bvec0)) 
# turn nvec into kxk matrix
  	 nmat <- matrix(0,nrow=k,ncol=k); count <- 1
	 for(i in 1:k){
      	   	 for (j in 1:k){
      	  	     if (i<=j){
      	     	     	nmat[i,j] <- nvec[count]
	     	     	count <- count+1}}}
#	 for(i in 1:k){
#      	  	 for (j in 1:k){
#      	  	     if (i>j){
#      	     	     	nmat[i,j] <- nmat[j,i]
#		     }
#		 }}
#cat(nmat,"\n")
	 mainsum <- rep(0,k)
	 for (i in 1:k){
	        mainsum[i] <- 2*nmat[i,i] + bvec0[i]
	 }
	 for (i in 1:k){
	    for (j in 1:k){
	    	if (i < j) mainsum[i] <- mainsum[i] + nmat[i,j]
		if (i > j) mainsum[i] <- mainsum[i] + nmat[j,i]
#		   ifelse (j>i,mainsum[i] <- mainsum[i] + nmat[i,j],mainsum[i] <- mainsum[i] + nmat[j,i])}
	    }
	 }
         nsumoff <- sum(nvec) - sum( diag(nmat))
	 term1 <- nsumoff*log(2) 
#cat("n: ",nvec,"nmat ",nmat,"nsumoff ",nsumoff,"term1 ",term1,"mainsum ",mainsum,"\n")
	 DirichNormHWE <- DirichNormHWE + term1 + sum(lgamma(mainsum)) - lgamma(sum(mainsum))
	 DirichNormHWE <- exp(DirichNormHWE)
	 DirichNormHWE
}


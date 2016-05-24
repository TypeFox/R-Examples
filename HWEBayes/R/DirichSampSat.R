DirichSampSat <-
function(nvec,bvec,nsim){
# posterior samples for the distinct probs:
  	   if (length(nvec) != length(bvec) ) stop("DirichSampSat: Dimension of nvec and bvec differ\n") 
    	   pvec <-  rdirichlet(nsim,nvec+bvec)
	   k <- .5*(-1+sqrt(1+8*length(nvec)))
	   pmat <- array(0,dim=c(nsim,k,k)); count <- 1
# turn these into a kxk matrix
	   for(i in 1:k){
      	   	 for (j in 1:k){
      	  	     if (i<=j){
      	     	     	pmat[,i,j] <- pvec[,count]
	     	     	count <- count+1}}}
	   count <- k+1
	   for(i in 1:k){
      	   	 for (j in 1:k){
      	  	     if (i>j){
      	     	     	pmat[,i,j] <- pmat[,j,i]
	     		count <- count+1}}}
# now the marginal probs:
	   pmarg <- matrix(0,nrow=nsim,ncol=k)
	   for (i in 1:k){
	       pmarg[,i] <- pmat[,i,i] 
	       for (j in 1:k){
	       	   if(i!=j){
			pmarg[,i] <- pmarg[,i] + .5*pmat[,i,j]}}}
# the fixation coeffs:
	   fixind <- array(0,dim=c(nsim,k,k))
	   for(i in 1:k){
	   	 for(j in 1:k){
		       if (i>j){
		       	  fixind[,i,j] <- 1-pmat[,i,j]/(2*pmarg[,i]*pmarg[,j])}}}
	   list(pvec=pvec,pmat=pmat,pmarg=pmarg,fixind=fixind)
}


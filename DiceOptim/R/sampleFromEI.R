ranperm <- function(X, N) order(runif(N))
randLHS <- function(n,k) {
        P <- matrix(nrow = n, ncol = k)
        P <- apply(P, 2, ranperm, N = n)
        P <- P - 1 + matrix(runif(n * k), nrow = n, ncol = k)
	return(P/n)
}

sampleFromEI <- function(model,minimization=TRUE,n= 1,initdistrib = NULL,lower=rep(0,model@d),upper=rep(1,model@d),T=NULL){
	
  #Generic function to sample from a density proportional to EI
  
  #### Arguments:
  # model : a km object
  # minimization : do we perform a global minimization or maximization ?
  # n : the total size of the sample distributed from a density proportional to EI
  # initdistrib : how these points are generated (default: "LHS" sequence)
  # d : dimension of the input set
  # lower, upper : arrays of size d with the lower (resp. upper) bounds in each dimension
  # T: the threshold used to compute the EI
  
  if(is.null(model)){
    print("Error in samplefromEI. A km object needs to be given")
    return(NULL)
  }
  
  d <- model@d
  
  if (length(lower) != length(upper) ){
    print("Error in samplefromEI : lower and upper must have the same length")
    return(NULL)
  }
  
  if(is.null(initdistrib)) {
	initial.points <- t(lower+t(randLHS(n=d*1000,k=d))*(upper-lower))
	} else {
	initial.points <- initdistrib
	}
	if(d==1) initial.points <- matrix(initial.points,ncol=1)
	#prediction on these initial candidate points
	predictions <- predict(object=model,newdata=initial.points,type="UK",checkNames=FALSE,cov.compute = FALSE,se.compute = TRUE,light.return = TRUE)
	  
  mn <- predictions$mean
  sn <- predictions$sd
	
  if(is.null(T)){
    if(minimization) T <- min(model@y)
    if(!minimization) T <- max(model@y)
  }
  
  uu <- (mn-T)/sn
  if(minimization) uu <- -1*uu
  ei <- sn*( uu*pnorm(uu) + dnorm(uu) )
  
  if ( minimization)  ei[is.nan(ei)] <- (T-mn)
  if (!minimization)  ei[is.nan(ei)] <- (mn-T)  
  ei <- ei*(ei>0)
  if(sum(ei!=0)<n) {
    maxEI <- max(ei)
    if (maxEI > 0) {
      ei <- ei/maxEI+1/(1000*d)
    } else {
      ei <- rep(1,(1000*d))
    }
  }
  my.indices <- sample(1:(1000*d), n, replace = FALSE, prob = ei)
  my.points <- initial.points[my.indices,]
		
	if(d==1) my.points <- matrix(my.points,ncol=1)
	if(n==1) my.points <- matrix(my.points,ncol=d)
    
	return(my.points)
}
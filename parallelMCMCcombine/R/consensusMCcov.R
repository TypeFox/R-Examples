consensusMCcov <- function( subchain, shuff = FALSE )
{
  
  ddata = length(dim(subchain))
  
  if (ddata !=3 )
  {
    stop("The subchain must be an array of dimension c(d,sampletotT,M).")    
  }
  
  d          <- dim(subchain)[1]  
  sampletotT <- dim(subchain)[2]    
  M          <- dim(subchain)[3]
  
  if ( M==1 )
  { 
    theta <- array(subchain[,,1],c(d,sampletotT))
    return (theta)
  }
  
  if (shuff == TRUE) # randomly permute subposteriors 
  {
    for (k in 1:M)
      subchain[,,k]<- subchain[,sample(sampletotT),k]
  }
  
  ## Assuming there is covariance for the model parameters
  # compute sigmahatm & sigmahatm.inverse (= W_s)
  
  sigmahatm     <- array(NA,c(d,d,M))    
  sigmahatM     <- matrix(NA,d,d)
  sigmahatM.pre <- matrix(NA,d,d)        
  
  if (d==1)
  {
    for (k in 1:M) { sigmahatm[1,,k] <- var(subchain[1,,k]) }    
  }
  else
  {
    for (k in 1:M) { sigmahatm[,,k] <- cov(t(subchain[,,k])) }
  }
  
  # computation of the inverses of covariance matrices (with try()):
  
  sigmahatm.inverse <- array(NA,dim=c(d,d,M))
    
  for (k in 1:M){
    
    res <- try( sigmahatm.inverse[,,k] <- solve(sigmahatm[,,k]), silent=TRUE)  
    
    if (class(res) == "try-error")
    {  
      stop(paste("Computation of the inverse of a covariance matrix for",
                 "one of the sample vectors in the data subset #",k,"is failed.",
                 "Here is the system R error-message:\n",attr(res,"condition")))  
    }
  }            

  sigmahatM <- solve(rowSums(sigmahatm.inverse, dims=2)) # (sum W_s)^(-1)
  
  
  # Compute unified posterior samples
  
  theta <- matrix(NA,nrow=d,ncol=sampletotT) # the resulting posterior samples

  wvec  <- array(NA, c(d,1))
  
  for (i in 1:sampletotT){

    wvec <- rep(0,d)
        
    for (s in 1:M){
      wvec <- wvec + sigmahatm.inverse[,,s] %*% subchain[,i,s] 
    }

    theta[,i] <- sigmahatM %*% wvec        
  }
  
  return(theta)
}

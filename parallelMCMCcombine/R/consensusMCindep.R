consensusMCindep <- function( subchain, shuff = FALSE )
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
  else
  {          
    if (shuff == TRUE) # randomly permute subposteriors 
    {
      for (s in 1:M)
        subchain[,,s]<- subchain[,sample(sampletotT),s]      
    }
    
    # compute weights W

    W <- matrix(NA,d,M)
    
    for (j in 1:d)
      for (s in 1:M)
        W[j,s] <- 1/var(subchain[j,,s])
    
    # compute unified posterior samples    
    
    theta <- matrix(NA,nrow = d,ncol=sampletotT)              
    
    for (i in 1:sampletotT)    
      theta[,i] <- rowSums( W * subchain[,i,], dims=1)/ rowSums(W)      
        
    return(theta)    
  }  
  
}

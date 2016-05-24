sampleAvg <- function(subchain, shuff=FALSE)
{
  # Check-1 Subposterior data must have the right format:
  # (subchain should be an array of size d x sampletotT x M)
  
  ddata = length(dim(subchain));
  
  if (ddata !=3 )
  {
    stop("The subchain must be an array of dimension c(d,sampletotT,M).")    
  }
  else
  {         
    d          <- dim(subchain)[1]  
    sampletotT <- dim(subchain)[2]    
    M          <- dim(subchain)[3]
    
    if (M==1) 
      { 
        theta <- array(subchain[,,1],c(d,sampletotT))        
        return (theta) 
      }    
    else
    {
      # shuffle subposterior sample vector for each data subset
      
      if (shuff == TRUE)
      {
        for (k in 1:M)
        {
          subchain[,,k]<- subchain[,sample(sampletotT),k]
        }
      }          
  
      if (d==1)
      {
        theta <- array( rowMeans(subchain[,,], dims=1), c(d,sampletotT))                        
      }
      else
      {
        theta <- array( rowMeans(subchain[,,], dims=2), c(d,sampletotT))
      }
            
      return (theta)
    }    
  }  
}
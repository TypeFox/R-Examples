
# Creates model matrix for quadratic multivariate regression 

.quadModMat <- function(X)
{
  nPar <- ncol(X)
  
  # Center X
  X <- t( t(X) - colMeans(X) )

  # Add quadratic terms
  mod <- cbind(1, X, (X^2)/2)
  
  # Add interactions
  if(nPar > 1){
    
    comb <- t( combn(nPar, 2) )
    
    for(jj in 1:nrow(comb)){ 
      
      mod <- cbind(mod, X[ , comb[jj, 1]] * X[ , comb[jj, 2]]) 
      
    }
    
  }
  
  return( mod ) 
} 
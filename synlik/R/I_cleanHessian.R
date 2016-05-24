# Cleans the hessian used to get confidence intervals for the parameters
# Until the hessian is not positive definite, we take the smallest eigenvalues
# and increase them to a tolerance. We then invert this modified hessian and 
# we remove the parameter with the highest variance from the hessian.
# We returned the reduced hessian and the indexes of the parameters that have
# been eliminated.

.cleanHessian <- function(hessian, tol = 1e-8)
{
  badParam <- numeric()
  
  repeat{
    
    eDec <- eigen(hessian)
    
    # Find small eigenvalues
    lowEigen <- which(eDec$values < eDec$values[1] * tol)
    
    # If there are none, I consider the hessian to be PD and I exit
    if( length(lowEigen) == 0 ) break
    
    # Increase the small eigenvalues to a threshold
    eDec$values[lowEigen] <- eDec$values[1] * tol
    
    # Invert the modified hessian to get the covariance
    COV <- eDec$vectors%*%(t(eDec$vectors) / eDec$values)
        
    # I identify the parameter with the highest variance, I remove the corresponding
    # elements from the Hessian and I store its index
    bad <- which.max( diag(COV) )
    hessian <- hessian[-bad, -bad]
    offset <- sum( badParam <= bad )
    
    badParam <- c(badParam, bad + offset)
    
  }
  
  return( list("hessian" = hessian, "badParam" = badParam) )
}
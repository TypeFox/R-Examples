##Transform correlation matrix to [0,1] space

transform.correlations <- function(correlations){
  #INPUT: correlations is a d x d matrix whose (i,j)th entry is the
  #       correlation between i and j
  #First the observed correlation matrix is transformed to partial correlations. Next,
  #this is transformed to the [0,1] space. 
  
  #Initialize
  d <- dim(correlations)[1]
  partials <- matrix(1, d, d)
  
  #NOTE: partials should be symmetric. If not, send an error
  if(isSymmetric(correlations) == FALSE){
    stop("The correlation matrix must be symmetric.")
  }
  
  if(sum(eigen(correlations)$values < 0) > 0){
    stop("The correlation matrix must be non-negative definite.")
  }
  #super-diagonal
  diag(partials[-nrow(partials), -1]) <- diag(correlations[-nrow(correlations), -1])
  
  #sub-diagonal
  diag(partials[-1, -ncol(partials)]) <- diag(correlations[-1, -ncol(correlations)])
  
  #Recursively obtain the remaining correlations using the transformation
  #given in the Harry Joe paper.
  
  for(k in 2 : (d - 1)){
    for(i in 1 : (d - k)){
      
      R2 <- correlations[(i + 1): (i + k - 1), (i + 1): (i + k - 1)]
      r1 <- correlations[i, (i + 1) : (i + k - 1)]
      r3 <- correlations[i + k, (i + 1): (i + k - 1)]
      
      D <- sqrt((1 - t(r1) %*% solve(R2) %*% r1) * (1 - t(r3) %*% solve(R2) %*% r3))
      
      partials[i, i + k] <- as.numeric((correlations[i, i + k] - t(r1) %*% solve(R2) %*% r3) / D)
      
      partials[i + k, i] <- partials[i, i + k]
    }
  }
  #transform the partial correlations to the [0,1] space for GERGM estimation
  transformed.data <- (partials + 1) / 2 
  return(transformed.data)
}



####Transform [0,1] space back to correlation space
bounded.to.correlations <- function(bounded.network){
  #transform back to the partial space
  partials <- 2*bounded.network - 1
  diag(partials) <- 1
  #transform to correlation space
  correlations <- partials.to.correlations(partials)
}

partials.to.correlations <- function(partials){
  #INPUT: partials is a d x d matrix whose (i,j)th entry is the
  #       sampled partial correlation between i and j given 
  #       all other terms in between i and j
  #NOTE:  partials must be a symmetric matrix. The super- and 
  #       sub-diagonal terms are exactly the one step correlations
  
  #Initialize
  d <- dim(partials)[1]
  correlations <- matrix(1, d, d)
  
  #NOTE: partials should be symmetric. If not, send an error
  if(isSymmetric(partials) == FALSE){
    stop("The partial correlation matrix must be symmetric.")
  }
  
  #super-diagonal
  diag(correlations[-nrow(correlations), -1]) <- diag(partials[-nrow(partials), -1])
  
  #sub-diagonal
  diag(correlations[-1, -ncol(correlations)]) <- diag(partials[-1, -ncol(partials)])
  
  #NOTE: the above should be symmetric. If not, send an error
  if(isSymmetric(correlations) == FALSE){
    stop("The correlation matrix must be symmetric. Check the super- and
         sub- diagonals of object partials.")
  }
  
  #Recursively obtain the remaining correlations using the transformation
  #given in the Harry Joe paper.
  
  for(k in 2 : (d - 1)){
    for(i in 1 : (d - k)){
      
      R2 <- correlations[(i + 1): (i + k - 1), (i + 1): (i + k - 1)]
      r1 <- correlations[i, (i + 1) : (i + k - 1)]
      r3 <- correlations[i + k, (i + 1): (i + k - 1)]
      
      D <- sqrt((1 - t(r1) %*% solve(R2) %*% r1)* (1 - t(r3) %*% solve(R2) %*% r3))
      correlations[i, i + k] <- t(r1) %*% solve(R2) %*% r3 + partials[i, i + k]*D 
      correlations[i + k, i] <- correlations[i, i + k]
    }
  }
  return(correlations)
}
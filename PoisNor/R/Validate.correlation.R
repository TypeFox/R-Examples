Validate.correlation <-
function(no.pois, no.norm, corMat, lamvec){
    
  nPois = length(lamvec)
  nNorm = ncol(corMat) - nPois
  
  samples=100000
  u = runif(samples, 0, 1)
  X = rnorm(samples,0,1)
  Y = rnorm(samples,0,1)
  U = pnorm(X)
  
  errorCount =0
  
  if ( ncol(corMat) != (no.pois + no.norm) ){
    stop("Dimension of correlation matrix does not match the number of variables!\n")
  } 
    
  if (is.positive.definite(corMat) == FALSE) {
    stop("Specified correlation matrix is not positive definite! \n")
  }
  if (isSymmetric(corMat) == FALSE) {
    stop("Specified correlation matrix is not symmetric! \n")
  }
  
  
  if ( sum(corMat>1)>0) {
    stop("Correlation values cannot be greater than 1! \n")
  }
  
  if ( sum(corMat<(-1) )>0) {
    stop("Correlation values cannot be less than -1! \n")
  }
  
  if ( sum(diag(corMat)!=1)>0) {
    stop("All diagonal elements of the correlation matrix must be 1! \n")
  }
  
  
  maxcor = Valid.correlation(no.pois, no.norm,lamvec)$max
  mincor = Valid.correlation(no.pois, no.norm,lamvec)$min
  
  for (i in 1:nrow(corMat)){
    for (j in 1:i){
      if (errorCount==0)  cat(".")
      if (i!=j){         
        if( corMat[i,j] >  maxcor[i,j] | corMat[i,j] < mincor[i,j] ) {
          cat("\n corMat[",i,",",j,"] must be between ", round(mincor[i,j],3), " and ",round(maxcor[i,j],3),"\n")  
          errorCount= errorCount + 1
        }
      }  
            
    }
      
  }
  
  cat("\n")
  if (errorCount>0) {  
   stop ("Range violation occurred in the target correlation matrix!\n")
  }
  
  return(TRUE)
}

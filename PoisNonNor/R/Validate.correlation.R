Validate.correlation <-
function (cmat, pmat=NULL, lamvec=NULL) 
{
  
  errorCount = 0
  if (isSymmetric(cmat) == FALSE) {
    stop("Specified correlation matrix is not symmetric! \n")
  }
  if (sum(cmat > 1) > 0) {
    stop("Correlation values cannot be greater than 1! \n")
  }
  if (sum(cmat < (-1)) > 0) {
    stop("Correlation values cannot be less than -1! \n")
  }
  if (sum(diag(cmat) != 1) > 0) {
    stop("All diagonal elements of the correlation matrix must be 1! \n")
  }
  if (!is.positive.definite(cmat)) {
    stop("correlation matrix must be PD \n")
  }
  
  n1 = ifelse(is.null(lamvec),0,length(lamvec)) 
  n2 = ifelse(is.null(pmat),0,dim(pmat)[1])
  
  if ((n1+n2) != dim(cmat)[1]) {
    stop("Correlation matrix dimension is not consistent with number of variables!\n")
  }
  
  if ( (!is.null(lamvec)) & (sum(lamvec>0) < n1) ) {
    stop("Specified lambda should be positive \n")
  }
  
  if ((!is.null(pmat)) & (dim(pmat)[2] != 4)){
    stop("column of pmat must be 4\n")
  }
  
  maxcor = diag(NA,(n1+n2))
  mincor = diag(NA,(n1+n2))
  
  if (n1 != 0) {
    bounds.corr.PP =   bounds.corr.GSC.PP (lamvec) 
    maxcor[1:n1,1:n1] = bounds.corr.PP$max
    mincor[1:n1,1:n1] = bounds.corr.PP$min 
    
  }  
  if (n2 != 0) {
    bounds.corr.NN = bounds.corr.GSC.NN (pmat)
    maxcor[(n1+1):(n1+n2),(n1+1):(n1+n2)] = bounds.corr.NN$max
    mincor[(n1+1):(n1+n2),(n1+1):(n1+n2)] = bounds.corr.NN$min
  }
    
  if (n1 != 0 & n2 != 0) {
    bounds.corr.NNP = bounds.corr.GSC.NNP (lamvec,pmat) 
    maxcor[1:n1,(n1+1):(n1+n2)] = bounds.corr.NNP$max
    maxcor[(n1+1):(n1+n2),1:n1] = t(maxcor[1:n1,(n1+1):(n1+n2)])
    mincor[1:n1,(n1+1):(n1+n2)] = bounds.corr.NNP$min
    mincor[(n1+1):(n1+n2),1:n1] = t(mincor[1:n1,(n1+1):(n1+n2)])
  }

  for (i in 2:(n1+n2)) {
    for (j in 1:(i-1)) {
      if (errorCount == 0)  cat(".")
      if (cmat[i, j] > maxcor[i, j] | cmat[i, j] < mincor[i, j]) {
        cat("\n cmat[", i, ",", j, "] must be between ", round(mincor[i, j], 3), " and ", round(maxcor[i, j], 3), "\n")
        errorCount = errorCount + 1
      }
    }
  }
  cat("\n")
  if (errorCount > 0) {
    stop("Range violation occurred in the target correlation matrix!\n")
  }
  return(TRUE)
}

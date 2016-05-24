RNG_P_NN <-
function(lamvec=NULL, cmat, rmat=NULL, norow, mean.vec = NULL, variance.vec = NULL){
  
  n1 = ifelse(is.null(lamvec),0,length(lamvec))
  n2 = ifelse(is.null(rmat),0,dim(rmat)[1]) 
  
  if ( (!is.null(lamvec)) & (sum(lamvec>0) < n1) ) {
    stop("Specified lambda should be positive \n")
  }
  
  if ((!is.null(rmat)) & (dim(rmat)[2] != 2)){
    stop("column of rmat must be 2\n")
  }
  
  if ( (!is.null(rmat)) & (sum(rmat[,2]>=(rmat[,1]^2-2)) < n2) ) {
    stop("Specified skewness and kurtosis parameter should be v2>=v1^2-2 \n")
  }
  
  if ((n1+n2) != dim(cmat)[1]) {
    stop("Correlation matrix dimension is not consistent with number of variables!\n")
  }
  
  cmat_N_N = diag(1,(n1+n2))
  
  pmat = NULL
  if (n2 != 0) {
    pmat = Param.fleishman(rmat)
  }


  if (Validate.correlation(cmat,pmat,lamvec)) cmat_N_N = intercor.all(cmat,pmat,lamvec)
  if (!is.positive.definite(cmat_N_N)) {
    warning("Intermediate correlation matrix is not positive definite. Nearest positive definite matrix is used!")
    cat(cmat_N_N,"\n")
    cmat_N_N = as.matrix(nearPD(cmat_N_N, corr = TRUE, keepDiag = TRUE)$mat)
    cat(cmat_N_N,"\n")
  }

  #RNG
  X = mvrnorm(n=norow,rep(0,dim(cmat)[1]),cmat_N_N)
  data = matrix(NA,nrow=norow,ncol=dim(cmat)[1])
  #transform poission part
  if (n1>0)  data[,1:n1]=t(qpois(t(pnorm(X[,1:n1])), lamvec))

  #transform NN part
  if (n2>0) {
    for (i in (n1+1):(n1+n2))
    {
      j = i - n1
      data[,i]= pmat[j,1]+pmat[j,2]*X[,i]+pmat[j,3]*X[,i]*X[,i]+pmat[j,4]*X[,i]*X[,i]*X[,i]
      if (!is.null(variance.vec)) data[,i] = data[,i]*sqrt(variance.vec[j])
      if (!is.null(mean.vec)) data[,i] = data[,i] + mean.vec[j]
    }
  }
  return(data)
}

intercor.all <-
function(cmat, pmat=NULL, lamvec=NULL){
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
  
  
  cmat_N_N = diag(1,(n1+n2))
  
  if (n1 != 0) {
    cmat_PP = cmat[1:n1,1:n1]
    cmat_N_N[1:n1,1:n1] = intercor.PP (lamvec,cmat_PP)
  }  
  if (n2 != 0) {
    cmat_NN = cmat[(n1+1):(n1+n2),(n1+1):(n1+n2)]
    cmat_N_N[(n1+1):(n1+n2),(n1+1):(n1+n2)] = intercor.NN(pmat,cmat_NN)
  }
  
  if (n1 != 0 & n2 != 0) {
    cmat_NNP = cmat[1:n1,(n1+1):(n1+n2)]
    cmat_N_N[1:n1,(n1+1):(n1+n2)] = intercor.NNP(lamvec,cmat_NNP,pmat)
    cmat_N_N[(n1+1):(n1+n2),1:n1] = t(cmat_N_N[1:n1,(n1+1):(n1+n2)])
  }
  
  return(round(cmat_N_N,3))
}

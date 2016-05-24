intercor.NNP <-
function(lamvec, cmat, pmat){
  if ((length(lamvec) != dim(cmat)[1])|(dim(pmat)[1] != dim(cmat)[2])) {
    stop("Correlation matrix dimension is not consistent with number of variables!\n")
  }
  
  if (sum(lamvec <= 0) > 0) {
    stop("lambda should be positive \n")
  }
  
  if (dim(pmat)[2] != 4){
    stop("column of pmat must be 4\n")
  }
  
  n1 = length(lamvec)
  n2 = dim(pmat)[1]
  norow = 1e+05
  
  cor_NN = matrix(NA,nrow=n1,ncol=n2)
  for (i in 1:n1){
    for (j in 1:n2){
      X = rnorm(norow, 0, 1)
      Y = rnorm(norow, 0, 1)
      U = pnorm(X)
      Xpois = qpois(U, lamvec[i])
      c = cor(Xpois[order(Xpois)], Y[order(Y)])/cor(X[order(X)],Y[order(Y)])
      cor_NN[i,j] = cmat[i,j]/c/(pmat[j,2]+3*pmat[j,4])
    }
  }
  return(round(cor_NN,3))
}

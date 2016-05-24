lambda0.matrix=function(x,lambda=0,maxit=100,trace.it=FALSE,thresh=1e-5){
  ismiss=is.na(x)
  if(any(ismiss))x[ismiss]=0
  svd(x)$d[1]
}

lambda0.Incomplete=function(x,lambda=0,maxit=100,trace.it=FALSE,thresh=1e-5){
  fit=svd.als(x,rank.max=2,lambda=lambda,maxit=maxit,trace.it=trace.it, thresh=thresh)
  lam0=fit$d[1]
  if(lam0<=thresh)warning("lambda too large; lambda0 is smaller than lambda")
  lam0+lambda
}
setGeneric("lambda0",lambda0.matrix)
setMethod("lambda0","Incomplete",lambda0.Incomplete)
setMethod("lambda0","SparseplusLowRank",lambda0.Incomplete)
setMethod("lambda0","sparseMatrix",lambda0.Incomplete)

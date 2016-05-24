Param.fleishman <-
function(rmat){
  if (dim(rmat)[2] != 2) {
    stop("column of rmat must be 2 \n")
  }
  
  if (sum(rmat[,2]>=(rmat[,1]^2-2)) < dim(rmat)[1]) {
    stop("Specified skewness and kurtosis parameter should be v2>=v1^2-2 \n")
  }
  
  pmat = matrix(NA,nrow=dim(rmat)[1],ncol=3)  
  for (i in 1:dim(rmat)[1]){
    pmat[i,] = BBsolve(par=rep(0,3), fn=fleishman.roots, r=rmat[i,])$par
  }
  pmat = as.matrix(cbind(-pmat[,2],pmat))
  
  return(round(pmat,7))
}

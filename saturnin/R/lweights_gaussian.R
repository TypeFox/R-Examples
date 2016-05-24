lweights_gaussian <-
function(data,a = ncol(data),mu = numeric(p),au = 1,T = diag(ncol(data),ncol(data)),nbcores=1){
  p <- ncol(data)
  n <- nrow(data)
  datamean <- apply(data,2,mean)
  R <- T + (n-1)*cov(data) + (au*n/(au+n))*matrix((mu - datamean),p,1)%*%matrix((mu - datamean),1,p)
  uptri <- upper.tri(matrix(0,p,p))
  
  if (requireNamespace("parallel",quietly = TRUE) &&
        (nbcores > 1)){
    lW <- matrix(parallel::mcmapply(function(y, i, j) if (y){
      0.5*(a-p+2)*ldet(T[c(i,j),c(i,j)]) - 0.5*(a-p+n+2)*ldet(R[c(i,j),c(i,j)])   
    } else 0,
    uptri, 
    row(uptri), 
    col(uptri),
    mc.cores = nbcores), 
    nrow = nrow(uptri))  
  } else {
    lW <- matrix(mapply(function(y, i, j) if (y){
      0.5*(a-p+2)*ldet(T[c(i,j),c(i,j)]) - 0.5*(a-p+n+2)*ldet(R[c(i,j),c(i,j)])
    } else 0,
    uptri, 
    row(uptri), 
    col(uptri)), 
    nrow = nrow(uptri))
  }

  diaglW <- sapply(1:p,function(i) 0.5*(a-p+1)*log(abs(T[i,i]))  -0.5*(a-p+n+1)*log(abs(R[i,i])))
  lW <- lW + t(lW)
  lW <- sapply(1:p,function(x) lW[,x] - diaglW[x])
  lW <- t(sapply(1:p,function(x) lW[x,] - diaglW[x]))
  diag(lW) <- 0
  return(lW)
}

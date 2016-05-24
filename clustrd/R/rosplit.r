rosplit <- function(data,U0)
{
  # cluster splitting algorithm 
  # n = number of objects
  # k = number of clusters of the partition
  # maxiter = max number of iterations
  out=list()
  maxiter = 99
  n = dim(U0)[1]
  k = dim(U0)[2]
  
  data = data.matrix(data)
  
  n = dim(data)[1]
  m = dim(data)[2]
  
  un = data.matrix(rep(1,n))
  uk = data.matrix(rep(1,k))
  um = data.matrix(rep(1,m))
  
  distt=matrix(0,k,n)
  eps=0.0000000001
  st=sum(sum(data^2))
  so=0
  Xuk = kronecker(data,uk)
  
  Xmean0 = pseudoinverse(U0)%*%data
  
  for (iter in c(1:maxiter))
  {
    # given Xmean0 assign each units to the closest cluster
    distt = matrix(((Xuk-kronecker(un,Xmean0))^2)%*%um,k,n)
    U = matrix(0,n,k)
    for (i in c(1:n)) {
      p = which.min(distt[,i])
      U[i,p] = 1
    }
    su = apply(U,2,sum)
    # given U compute Xmean (compute centroids)
    Xmean = pseudoinverse(U)%*%data
    
    # stopping rule
    sa = 100*sum(sum((U%*%Xmean)^2))/st
    dif=sa-so
    
    if ((dif > eps) && (sum(sum(abs(U-U0))) > 1)) {
      Xmean0=Xmean
      U0=U
      so=sa
    } else {
      break
    }
  }
  out$U = U
  out
}
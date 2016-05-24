transform_z <- function(data,is.weight=TRUE,is.exact=TRUE,r=1,c=1){
  
  ## function that transforms a categorical data matrix in a standardized residual matrix accroding to 
  ## the weights in D1 e D2
 # require(dummies)
  out=list()
  
  data=data.frame(data)
  data=data.frame(lapply(data,as.factor))
  ddZ=as.matrix(dummy.data.frame((data),drop=F))

  out$dZ=ddZ
  out$J=ncol(out$dZ)
  out$Q=ncol(data)
  
  if(is.weight==T && is.exact==F) {
    n=nrow(data)
    dZ=out$dZ  
    CZ=dZ/sum(dZ)
    r=apply(CZ,1,sum)
    eCZ1=r[1:n] %*% t(c)
    c=apply(CZ,2,sum)
    eCZ2=r[1:n] %*% t(c)
    SZ = (CZ-eCZ1)/sqrt(eCZ2)
    out$r=r
    out$c=c
  }
  
  if(is.weight==T && is.exact==T) {
    n=length(r)
    n1=nrow(data)
    Q=ncol(data)
    dZ=out$dZ  
    CZ=dZ/(n*Q)
    eCZ=r[1:n1] %*% t(c)
    SZ = (CZ-eCZ)/sqrt(eCZ)
    out$r=r
    out$c=c
  }
  
  
  
  if(is.weight==F) {
    n=nrow(data)
    Q=ncol(data)
    dZ=out$dZ  
    CZ=dZ/sum(dZ)
    r=apply(CZ,1,sum)
    c=apply(CZ,2,sum)
    eCZ=r[1:n] %*% t(c)
    SZ = (CZ-eCZ)/sqrt(eCZ)
    out$r=r
    out$c=c
  }
  
  
  out$SZ=SZ
  out
}
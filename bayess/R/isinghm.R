isinghm=function(niter,n,m=n,beta){
  x=sample(c(0,1),n*m,prob=c(0.5,0.5),replace=TRUE)
  x=matrix(x,n,m)
  for (i in 1:niter){
    sampl1=sample(1:n)
    sampl2=sample(1:m)
    for (k in 1:n){
    for (l in 1:m){
     n0=xneig4(x,sampl1[k],sampl2[l],x[sampl1[k],sampl2[l]])
     n1=xneig4(x,sampl1[k],sampl2[l],1-x[sampl1[k],sampl2[l]])
     if (runif(1)<exp(beta*(n1-n0)))
       x[sampl1[k],sampl2[l]]=1-x[sampl1[k],sampl2[l]]
    }}}
  x
  }

getstartvals<-function(x,eps=1e-03)
  {
    # Returns starting values for the relative precision matrix delta 
    n<-ncol(x)
    startvals<-double(n+n*(n+1)/2)
    startvals[1:n]<-apply(x,2,mean,na.rm=TRUE)
    
    sampmat<-cov(x,use="p") # sample var-cov matrix
    eig<-eigen(sampmat,symmetric=TRUE)
    realvals<-sapply(eig$values, function(y) ifelse(is.complex(y),0,y))
    smalleval<-eps*min(realvals[realvals>0])
    posvals<-pmax(smalleval,realvals)
    mypdmat<-eig$vectors %*% diag(posvals) %*% t(eig$vectors)
    myfact<-chol(mypdmat)
    mydel<-solve(myfact,diag(n))
    signchange<-diag(ifelse(diag(mydel)>0,1,-1))
    mydel<-mydel %*% signchange # ensure that diagonal elts are positive                
    startvals[(n+1):(2*n)]<-log(diag(mydel))
    for(i in 2:n){   # assume n>2
      startvals[(2*n+sum(1:(i-1))-i+2):(2*n+sum(1:(i-1)))]<-mydel[1:(i-1),i]
    }
    startvals
  }

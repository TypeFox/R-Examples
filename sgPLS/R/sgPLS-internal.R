normv <- function(x) sqrt(sum(x**2))

soft.thresholding <- function(x,lambda){
  tol <- .Machine$double.eps ^ 0.5 
  y <- abs(x)-lambda 
  test  <- y < tol
  return(sign(x)*y*(1-test))  
}  

 

soft.thresholding.group <- function(x,ind,lambda){
  tab.ind <- c(0,ind,length(x))
  tol <- .Machine$double.eps ^ 0.5 
  res <- NULL
  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]  
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    y <- 1-(lambda/2)*sqrt(ji)/normv(vecx)
    if(y < tol) y <- 0
    res <- c(res,vecx*y)  
  }
  return(res)    
}



soft.thresholding.sparse.group <- function(x,ind,lambda,alpha,ind.block.zero){
  tab.ind <- c(0,ind,length(x))
  res <- NULL
  for (i in 1:(length(ind)+1)){
    ji <- tab.ind[i+1]-tab.ind[i]  
    vecx <- x[((tab.ind[i]+1):tab.ind[i+1])]
    if(i%in%ind.block.zero) {vecx <- rep(0,ji)} else{
      temp <- soft.thresholding(vecx,lambda*alpha/2)
      vecx <- 0.5*temp*(1-lambda*(1-alpha)*sqrt(length(vecx))/sqrt(sum(temp**2))) 
    }
    res <- c(res,vecx)  
  }
  return(res)    
}



lambda.quadra <- function(x,vec,alpha){
  return(sum(soft.thresholding(vec,x*alpha/2)**2)-length(vec)*((1-alpha)*x)**2)
}

step1.spls.sparsity <- function(X,Y,sparsity.x,sparsity.y,epsilon,iter.max){
  n <- dim(X)[1]
  Z <- t(X)%*%Y
  svd.Z <- svd(Z,nu=1,nv=1)
  
  u.tild.old <- svd.Z$u
  v.tild.old <- svd.Z$v
  u.tild.previous <- v.tild.previous <- 0
  iter <- 0
  #|(norm(v.tild.old-v.tild.previous)>epsilon)) TO BE ADDED PROBLEM IN MIXOMICS ???
  ### Step c
  while (((normv(u.tild.old-u.tild.previous)/normv(u.tild.old))>epsilon)  & (iter <iter.max))  {
    if(sparsity.x==0) {lambda.x <- 0} else{
      lambda.x <- sort(abs(Z%*%matrix(v.tild.old,ncol=1)))[sparsity.x]}
    u.tild.new <- soft.thresholding(Z%*%matrix(v.tild.old,ncol=1),lambda=lambda.x)
    u.tild.new <- u.tild.new/sqrt(sum(u.tild.new**2))
    if(sparsity.y==0) lambda.y <- 0 else lambda.y <- sort(abs(t(Z)%*%matrix(u.tild.old,ncol=1)))[sparsity.y]
    v.tild.new <- soft.thresholding(t(Z)%*%matrix(u.tild.old,ncol=1),lambda=lambda.y)
    v.tild.new <- v.tild.new/sqrt(sum(v.tild.new**2))
    
    u.tild.previous <- u.tild.old
    v.tild.previous <- v.tild.old
    
    u.tild.old <- u.tild.new
    v.tild.old <- v.tild.new
    
    iter <- iter +1
  }  
  res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 
  
}


step1.sparse.group.spls.sparsity <- function(X,Y,ind.block.x,ind.block.y,sparsity.x,sparsity.y,epsilon,iter.max,alpha.x,alpha.y,upper.lambda=upper.lambda){
  n <- dim(X)[1]
  Z <- t(X)%*%Y
  svd.Z <- svd(Z,nu=1,nv=1)
  
  u.tild.old <- svd.Z$u
  v.tild.old <- svd.Z$v
  u.tild.previous <- v.tild.previous <- 0
  iter <- 0
  
  ### Step c
  #|(norm(v.tild.old-v.tild.previous)>epsilon)
  while (((normv(u.tild.old-u.tild.previous)>epsilon) ) & (iter <iter.max))  {
    vecZV <- Z%*%matrix(v.tild.old,ncol=1)
    tab.ind <- c(0,ind.block.x,length(vecZV))
    lamb.x <- NULL
    lamb.max <- upper.lambda
    for (i in 1:(length(ind.block.x)+1)){
      ji <- tab.ind[i+1]-tab.ind[i]  
      vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
      lamb.x <- c(lamb.x,uniroot(lambda.quadra,lower=0,upper=lamb.max,vec=vecx,alpha=alpha.x)$root)
    }   
    if(sparsity.x==0){lambda.x <- sort(lamb.x)[1]-1} else {
      lambda.x <- sort(lamb.x)[sparsity.x]}
    
    ####block to zero
    index.block.zero.x <- which(lamb.x<=lambda.x)
    
    
    if(sparsity.y==0) {lambda.y <- 0} else { 
      vecZV <- t(Z)%*%matrix(u.tild.old,ncol=1)
      tab.ind <- c(0,ind.block.y,length(vecZV))
      lamb.y <- NULL
      lamb.max <- 100000
      res <- NULL
      for (i in 1:(length(ind.block.y)+1)){
        ji <- tab.ind[i+1]-tab.ind[i]  
        vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
        lamb.y <- c(lamb.y,uniroot(lambda.quadra,lower=0,upper=lamb.max,vec=vecx,alpha=alpha.y)$root)
      }
      lambda.y <- sort(lamb.y)[sparsity.y]
      index.block.zero.y <- which(lamb.y<=lambda.y)
    }
    
    u.tild.new <- soft.thresholding.sparse.group(Z%*%matrix(v.tild.old,ncol=1),ind=ind.block.x,lambda=lambda.x,alpha=alpha.x,ind.block.zero=index.block.zero.x)
    
    u.tild.new <- u.tild.new/sqrt(sum(u.tild.new**2))
    
    if(sparsity.y==0) {v.tild.new <- t(Z)%*%matrix(u.tild.old,ncol=1)} else {
      v.tild.new <- soft.thresholding.sparse.group(t(Z)%*%matrix(u.tild.old,ncol=1),ind=ind.block.y,lambda=lambda.y,alpha=alpha.y,ind.block.zero=index.block.zero.y)
    }
    
    v.tild.new <- v.tild.new/sqrt(sum(v.tild.new**2))
  
    u.tild.previous <- u.tild.old
    v.tild.previous <- v.tild.old
    
    u.tild.old <- u.tild.new
    v.tild.old <- v.tild.new
    
    iter <- iter +1
  }  
  res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 
  
}




step1.group.spls.sparsity <- function(X,Y,ind.block.x,ind.block.y,sparsity.x,sparsity.y,epsilon,iter.max){
  n <- dim(X)[1]
  Z <- t(X)%*%Y
  svd.Z <- svd(Z,nu=1,nv=1)
  
  u.tild.old <- svd.Z$u
  v.tild.old <- svd.Z$v
  u.tild.previous <- v.tild.previous <- 0
  iter <- 0
  
  ### Step c
  #|(norm(v.tild.old-v.tild.previous)>epsilon)
  while (((normv(u.tild.old-u.tild.previous)>epsilon) ) & (iter <iter.max))  {
    
    vecZV <- Z%*%matrix(v.tild.old,ncol=1)
    tab.ind <- c(0,ind.block.x,length(vecZV))
    res <- NULL
    for (i in 1:(length(ind.block.x)+1)){
      ji <- tab.ind[i+1]-tab.ind[i]  
      vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
      res <- c(res,2*normv(vecx)/sqrt(ji))
    }
    if(sparsity.x==0) lambda.x <- 0 else{
      lambda.x <- sort(res)[sparsity.x]}
    
    
    if(sparsity.y==0) {lambda.y <- 0} else { 
      vecZV <- t(Z)%*%matrix(u.tild.old,ncol=1)
      tab.ind <- c(0,ind.block.y,length(vecZV))
      res <- NULL
      for (i in 1:(length(ind.block.y)+1)){
        ji <- tab.ind[i+1]-tab.ind[i]  
        vecx <- vecZV[((tab.ind[i]+1):tab.ind[i+1])]
        res <- c(res,2*normv(vecx)/sqrt(ji))
      }
      
      lambda.y <- sort(res)[sparsity.y]}
    
    
    u.tild.new <- soft.thresholding.group(Z%*%matrix(v.tild.old,ncol=1),ind=ind.block.x,lambda=lambda.x)
    u.tild.new <- u.tild.new/sqrt(sum(u.tild.new**2))
    
    if(sparsity.y==0) {v.tild.new <- t(Z)%*%matrix(u.tild.old,ncol=1)} else {
      v.tild.new <- soft.thresholding.group(t(Z)%*%matrix(u.tild.old,ncol=1),ind=ind.block.y,lambda=lambda.y)}
    
    v.tild.new <- v.tild.new/sqrt(sum(v.tild.new**2))
       
    u.tild.previous <- u.tild.old
    v.tild.previous <- v.tild.old
    
    u.tild.old <- u.tild.new
    v.tild.old <- v.tild.new
    
    iter <- iter +1
  }  
  res <- list(iter=iter, u.tild.new=u.tild.new,v.tild.new=v.tild.new) 
  
}


step2.spls <- function(X,Y,u.tild.new,v.tild.new,mode){
  ### Step d
  xi.h <- X%*% matrix(u.tild.new,ncol=1)/((normv(u.tild.new))**2)
  w.h  <- Y%*% matrix(v.tild.new,ncol=1)/((normv(v.tild.new))**2)
  
  ### Step e
  c.h <- t(X)%*%matrix(xi.h,ncol=1)/((normv(xi.h))**2)
  
  d.rh <- t(Y)%*%matrix(xi.h,ncol=1)/(sum(xi.h*xi.h))
  
  d.h <- t(Y)%*%matrix(w.h,ncol=1)/(sum(w.h*w.h))
  
  ###Step f and g
  X.h <- X - xi.h%*%t(c.h)
  if (mode=="regression") Y.h <- Y - xi.h%*%t(d.rh) else Y.h <- Y - w.h%*%t(d.h)
  
  res <- list(X.h=X.h,Y.h=Y.h,c=c.h,d=d.rh,e=d.h)
  return(res)
}



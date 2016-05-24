gPLS <- function(X,Y,ncomp,mode="regression",max.iter=500,tol=1e-06,keepX ,keepY=NULL,ind.block.x,ind.block.y=NULL){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  q <- ncol(Y)
  p <- ncol(X)
  n <- nrow(X)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    Y.names = "Y"
  else {
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names)) 
      Y.names = paste("Y", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Y)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }
  
  mat.c <-matrix(nrow = p, ncol = ncomp)
  mat.d <- matrix(nrow = q, ncol = ncomp)
  mat.e <- matrix(nrow = q, ncol = ncomp)
  
  
  X.s <- scale(X)
  Y.s <- scale(Y)
  
  sparsity.x <- length(ind.block.x)+1-keepX
  if(is.null(ind.block.y)) {sparsity.y <- rep(0,ncomp)} else {
    if (is.null(keepY)) keepY <- rep(length(ind.block.y)+1,ncomp)
    sparsity.y <- length(ind.block.y)+1-keepY}
  
  
  res.load <- step1.group.spls.sparsity(X=X.s,Y.s,ind.block.x=ind.block.x,ind.block.y=ind.block.y,sparsity.x=sparsity.x[1],sparsity.y=sparsity.y[1],epsilon=tol,iter.max=max.iter)
  res.deflat <- step2.spls(X=X.s,Y=Y.s,res.load$u.tild.new,res.load$v.tild.new,mode=mode)
  
  mat.c[,1] <- res.deflat$c
  iter <- res.load$iter
  if (mode=="regression") mat.d[,1] <- res.deflat$d else mat.e[,1] <- res.deflat$e
  load.u <- res.load$u.tild.new
  load.v <- res.load$v.tild.new
  
  if(ncomp>1) {
    
    for (h in 2:ncomp) {
      res.load <- step1.group.spls.sparsity(X=res.deflat$X.h,Y=res.deflat$Y.h,ind.block.x=ind.block.x,ind.block.y=ind.block.y,sparsity.x=sparsity.x[h],sparsity.y=sparsity.y[h],epsilon=tol,iter.max=max.iter)
      load.u <- cbind(load.u,res.load$u.tild.new)
      load.v <- cbind(load.v,res.load$v.tild.new)
      res.deflat <- step2.spls(X=res.deflat$X.h,Y=res.deflat$Y.h,res.load$u.tild.new,res.load$v.tild.new,mode=mode)
      mat.c[,h] <- res.deflat$c
      if (mode=="regression") mat.d[,h] <- res.deflat$d else mat.e[,h] <- res.deflat$e
      iter <- c(iter,res.load$iter)}
  }else{
    load.u <- matrix(load.u,ncol=1)
    load.v <- matrix(load.v,ncol=1)
  }
  
  mat.t <- X.s%*%load.u
  mat.u <- Y.s%*%load.v
       # mat.c regressor for X
       # mat.d regressor for Y "regression mode"
       # mat.e regressor for Y "canonical mode"
  cl = match.call()
  if (is.null(keepY)){
    if (is.null(ind.block.y)) keepY <- rep(ncol(Y),ncomp) else keepY <- rep(length(ind.block.y)+1,ncomp)
  } 
  
  result <- list(call = cl,X=X.s,Y=Y.s,ncomp=ncomp,mode=mode,keepX=keepX,keepY=keepY,mat.c=mat.c,mat.d=mat.d,mat.e=mat.e,loadings = list(X = load.u, Y = load.v),variates = list(X = mat.t, Y = mat.u),
                 names = list(X = X.names,Y = Y.names, indiv = ind.names),tol=tol,max.iter=max.iter,iter=iter,ind.block.x=ind.block.x,ind.block.y=ind.block.y)
  class(result) = c("gPLS","sPLS","spls", "pls")
  return(invisible(result))
}
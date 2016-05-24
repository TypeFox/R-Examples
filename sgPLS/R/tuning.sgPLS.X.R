tuning.sgPLS.X <- function(X,Y,folds=10,validation=c("Mfold","loo"),ncomp,keepX=NULL,alpha.x=NULL,grid.gX,grid.alpha.X,setseed,progressBar=FALSE,ind.block.x=ind.block.x,upper.lambda=10^9){
  choicesetseed <- setseed
  if(length(keepX)>(ncomp-1)) stop("The length of keepX should be less then ncomp")
  k <- 0
  l <- 0
  res <- matrix(0,ncol=length(grid.alpha.X),nrow=length(grid.gX))
  for (i in grid.gX){
    k <- k+1
    for (j in grid.alpha.X){
      if(is.null(keepX)) keepX1 <- rep(i,ncomp) else keepX1 <- c(keepX,rep(i,ncomp-length(keepX)))
      if(is.null(alpha.x)) alpha.x1 <- rep(j,ncomp) else alpha.x1 <- c(alpha.x,rep(j,ncomp-length(alpha.x)))
      l <- l+1
      cond <- TRUE
      while (cond) {
        model.sgpls <- sgPLS(X,Y,ncomp=ncomp,mode="regression",keepX=keepX1,ind.block.x=ind.block.x,alpha.x=alpha.x1,upper.lambda=upper.lambda)
        res.perf.sgpls <- try(perf(model.sgpls,criterion="MSEP",validation=validation,folds = folds,setseed=choicesetseed,progressBar=progressBar),silent=FALSE)
        if (class(res.perf.sgpls)[1]=="try-error"){ cond <- TRUE;choicesetseed=choicesetseed+1 } else {cond <- FALSE}
      }  
      res[k,l] <- sum(res.perf.sgpls$MSEP[,ncomp])
    }
    l <- 0
  }
  ind <- which.min(res)
  ind.XY <- which(res==res[ind],arr.ind=TRUE)
  keepX <- grid.gX[ind.XY[1]]
  alphaX <- grid.alpha.X[ind.XY[2]]
  return(list(MSEP=res,keepX=keepX,alphaX=alphaX))
}


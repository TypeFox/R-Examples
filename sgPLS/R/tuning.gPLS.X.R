tuning.gPLS.X <- function(X,Y,folds=10,validation=c("Mfold","loo"),ncomp,keepX=NULL,grid.X,setseed,progressBar=FALSE,ind.block.x=ind.block.x){
  choicesetseed <- setseed
  if(length(keepX)>(ncomp-1)) stop("The length of keepX should be less then ncomp")
  k <- 0
  res <- rep(0,length(grid.X))
  for (i in grid.X){
    if(is.null(keepX)) keepX1 <- rep(i,ncomp) else keepX1 <- c(keepX,rep(i,ncomp-length(keepX)))
    k <- k+1
    cond <- TRUE
    while (cond) {
      model.gpls <- gPLS(X,Y,ncomp=ncomp,mode="regression",keepX=keepX1,ind.block.x=ind.block.x)
      res.perf.gpls <- try(perf(model.gpls,criterion="MSEP",validation=validation,folds = folds,setseed=choicesetseed,progressBar=progressBar),silent=FALSE)
      if (class(res.perf.gpls)[1]=="try-error"){ cond <- TRUE;choicesetseed=choicesetseed+1 } else {cond <- FALSE}
    }  
    res[k] <- sum(res.perf.gpls$MSEP[,ncomp])
  }
  
  ind <- which.min(res)
  keepX <- grid.X[ind]
  return(list(MSEP=res,keepX=keepX))
}
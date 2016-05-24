tuning.sPLS.X <- function(X,Y,folds=10,validation=c("Mfold","loo"),ncomp,keepX=NULL,grid.X,setseed,progressBar=FALSE){
  choicesetseed <- setseed
  if(length(keepX)>(ncomp-1)) stop("The length of keepX should be less then ncomp")
  k <- 0
  res <- rep(0,length(grid.X))
  for (i in grid.X){
    if(is.null(keepX)) keepX1 <- rep(i,ncomp) else keepX1 <- c(keepX,rep(i,ncomp-length(keepX)))
    k <- k+1
    cond <- TRUE
    while (cond) {
      model.spls <- sPLS(X,Y,ncomp=ncomp,mode="regression",keepX=keepX1)
      res.perf.spls <- try(perf(model.spls,criterion="MSEP",validation=validation,folds = folds,setseed=choicesetseed,progressBar=progressBar),silent=FALSE)
      if (class(res.perf.spls)[1]=="try-error"){ cond <- TRUE;choicesetseed=choicesetseed+1 } else {cond <- FALSE}
    }  
    res[k] <- sum(res.perf.spls$MSEP[,ncomp])
  }
  
  ind <- which.min(res)
  keepX <- grid.X[ind]
  return(list(MSEP=res,keepX=keepX))
}
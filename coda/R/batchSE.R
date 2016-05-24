
"batchSE" <-
  function(x,batchSize=100) {
  UseMethod("batchSE")
}

"batchSE.mcmc" <- function(x,batchSize=100) {
  niter <- niter(x)
  nbatch <- niter%/%batchSize
  ## Truncate the odd lot observations
  ## Do this off the front instead of the back??
  niter <- nbatch*batchSize
  ibatch <- rep(1:nbatch,each=batchSize)[1:niter]
  batchMeans <- t(sapply(split(data.frame(x[1:niter,]),ibatch),
                               function(batch) apply(batch,2,mean)))
  grandMean <- apply(batchMeans,2,mean)
  mi2 <- sweep(batchMeans,2,grandMean,"-")^2
  stds<-sqrt(apply(mi2,2,sum)*batchSize/(nbatch-1))
  names(stds) <- dimnames(x)[[2]]
  stds/sqrt(niter(x))
}

"batchSE.mcmc.list" <- function(x,batchSize=100) {
  nchain <- nchain(x)
  niter <- niter(x)
  nbatch <- niter%/%batchSize
  ## truncate odd lot observations
  niter <- nbatch*batchSize
  ibatch <- rep(1:nbatch,each=batchSize)[1:niter]
  batchMeans <- NULL
  for (i in 1:nchain) {
    batchMeans <- rbind(batchMeans,
                        t(sapply(split(data.frame(x[[i]][1:niter,]),ibatch),
                                 function(batch) apply(batch,2,mean))))
  }
  #print(batchMeans)
  grandMean <- apply(batchMeans,2,mean)
  #cat("Grand Mean = ",grandMean,"\n")
  mi2 <- sweep(batchMeans,2,grandMean,"-")^2
  stds<-sqrt(apply(mi2,2,sum)*batchSize/(nchain*nbatch-1))
  names(stds) <- dimnames(x[[1]])[[2]]
  stds/sqrt(niter(x)*nchain(x))
}

## Needed for this function, but generally useful anyway.
as.data.frame.mcmc <- function(x, row.names = NULL, optional=FALSE, ...) {
    if (is.matrix(x))
        as.data.frame.matrix(x,row.names,optional, ...)
    else {
        if (is.null(row.names))
            row.names <- time(x)
        data.frame("var1"=as.numeric(x), row.names=row.names)
    }

}

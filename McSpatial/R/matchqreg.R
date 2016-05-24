matchqreg <- function(form,taumat=c(.10,.25,.50,.75,.90), qreglwr.smooth=TRUE, window=.50,bandwidth=0,kern="tcub", alldata=FALSE,
  graph.yhat=TRUE,graph.mean=TRUE,data) {

  mat <- model.frame(form,data=data)
  n = length(mat[,1])
  ntau <- length(taumat)
  timevect <- levels(factor(mat[,2]))
  ntime = length(timevect)
  yname <- names(mat)[1]
  xname <- names(mat)[2]

  if (qreglwr.smooth==FALSE) {
    yhat <- array(0,dim=c(ntime,ntau))
    for (i in seq(1:ntime)) {
      yhat[i,] <- quantile(mat[mat[,2]==timevect[i],1], taumat)
    }
  }
  if (qreglwr.smooth==TRUE) {
    fit <- qreglwr(mat[,1]~mat[,2],window=window,bandwidth=bandwidth,kern=kern,taumat=taumat)
    yhat <- aggregate(fit$yhat,by=list(mat[,2]),mean)
    yhat <- yhat[,-1]
  }
  rownames(yhat) <- timevect
  colnames(yhat) <- paste("q", gsub("0.","",taumat),sep="")
  
  if (graph.yhat==TRUE) {
    ymin = min(yhat)
    ymax = max(yhat)
    plot(timevect,yhat[,1],type="l",xlab=xname,ylab=yname,ylim=c(ymin,ymax),main="Quantiles")
    if (ntau>1) {
      for (j in seq(2,ntau)) {
        lines(timevect,yhat[,j])
      }
    }
  }
 
  ymean <- array(0,dim=ntime)
  for (i in seq(1:ntime)) {
    ymean[i] = mean(mat[mat[,2]==timevect[i],1])
  }
  if (graph.mean==TRUE) { plot(timevect,ymean,type="l",xlab=xname,ylab=yname,main="Means") }

  out <- list(yhat,ymean,timevect)
  names(out) <- c("yhat","ymean","timevect")
  return(out)
}


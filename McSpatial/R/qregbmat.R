qregbmat <- function(form,taumat=seq(.10,.90,.10),graphb=TRUE,graph.factor=FALSE,data=NULL) {

  xname <- colnames(model.matrix(form,data=data))
  ntau = length(taumat)
  nk = length(xname)
  bmat <- array(0,dim=c(ntau,nk))

  for (i in seq(1:ntau)) {
    fit <- rq(form,data=data,tau=taumat[i])
    bmat[i,] <- fit$coef
  }
    
  if (graphb==TRUE) {
    for (j in seq(1:nk)) {
      if (!(graph.factor==FALSE&substr(xname[j],1,6)=="factor" )) {plot(taumat,bmat[,j],xlab="Quantile",
        ylab="Coefficients",main=xname[j],type="l")}
    }
  }

  xname <- colnames(model.matrix(form,data=data))
#  xname[1] = colnames(model.matrix(form,data=data))[1]
  colnames(bmat) <- xname
  rownames(bmat) <- taumat
  return(bmat)
}



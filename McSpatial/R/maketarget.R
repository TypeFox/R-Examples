maketarget <- function(form,window=.25,bandwidth=0,kern="tcub",actualobs=FALSE,data=NULL){
  xmat <- as.matrix(model.matrix(form,data=data)[,-1])
  nk = ncol(xmat)

  if (nk==1&window>0)    {fit <- locfit(~lp(xmat[,1],nn=window),kern=kern) }
  if (nk==2&window>0)    {fit <- locfit(~lp(xmat[,1],xmat[,2],nn=window),kern=kern) }
  if (nk==1&bandwidth>0) {fit <- locfit(~lp(xmat[,1],h=2*bandwidth),kern=kern) }
  if (nk==2&bandwidth>0) {fit <- locfit(~lp(xmat[,1],xmat[,2],h=2*bandwidth),kern=kern) }
  xev <- lfeval(fit)$xev
  nt = length(xev)/nk
  target <- t(array(xev,dim=c(nk,nt)))
  target <- as.matrix(target)
  obs <- NULL

  if (actualobs==TRUE){
    if (nk==1) {vxmat <- var(xmat) }
    if (nk==2) {
      vxmat <- cov(xmat) 
    }
    obs <- array(0,dim=nt)
    for (i in seq(1:nt)) {
      dist <- sqrt(mahalanobis(xmat, target[i,], vxmat))
      obs[i] <- which.min(dist)
    }
    obs <- sort(unique(c(obs,chull(xmat))))
    target <- as.matrix(xmat[obs,])
  }

  nt = nrow(target)
  colnames(target) <- colnames(xmat)
  
  out <- list(target,obs)
  names(out) <- c("target","obs")
  return(out)
}


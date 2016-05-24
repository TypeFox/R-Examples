geogravity <- function(x,longvar,latvar,alpha=1,maxd=NULL,alldata=FALSE,window=.10,outmatrix=FALSE) {

  latvar  <- 2*pi*latvar/360
  longvar <- 2*pi*longvar/360

  if (alldata==TRUE) {targetobs <- seq(1,length(x)) }
  if (alldata==FALSE) {
    fit <- locfit(~lp(longvar,latvar,nn=window,deg=1)) 
    zev <- lfeval(fit)$xev
    nt = length(zev)/2
    target <- t(array(zev,dim=c(2,nt)))
    obs <- array(0,dim=nt)
    for (i in seq(1:nt)) {
      dist <- geodistance(longvar,latvar,target[i,1],target[i,2])$dist 
      obs[i] <- which.min(dist)
    }
    targetobs <- sort(unique(c(obs,chull(cbind(longvar,latvar)))))
  }
  nt = length(targetobs)
  
  gtarget <- array(0,dim=nt)
  dmat <- NULL
  if (outmatrix==TRUE&alldata==TRUE){dmat <- array(0,dim=nt,nt)}

  for (i in seq(1,nt)) {
    ii = targetobs[i]
    dist <- geodistance(longvar,latvar,longvar[ii],latvar[ii])$dist
    if (identical(maxd,NULL)) {dist <- x[ii]*x/(dist^alpha)}
    if (!identical(maxd,NULL)) {dist <- ifelse(dist<=maxd,  x[ii]*x/(dist^alpha), 0)}
    gtarget[i] = mean(dist[-ii],na.rm=TRUE)
    if (outmatrix==TRUE&alldata==TRUE){dmat[i,] <- dist}
  }
  if (outmatrix==TRUE&alldata==TRUE){diag(dmat) <- 0}

  gravity <- gtarget
  if (alldata==FALSE) {
    lmat <- cbind(longvar,latvar)
    gravity <- smooth12(lmat[targetobs,],gtarget, lmat)
  }
    
  out <- list(targetobs,gravity,dmat)
  names(out) <- c("targetobs","gravity","dmat")
  return(out)
}



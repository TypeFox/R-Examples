computeN <- function(x,lev){
  dat.cens <- matrix(x[x[,2]==1,],ncol=2)
  dat.uncens <- matrix(x[x[,2]==2,],ncol=2)
  
  n.cens <- as.vector(table(factor(dat.cens[,1][[1]],levels=lev)))
  n.uncens <- as.vector(table(factor(dat.uncens[,1][[1]],levels=lev)))
  
  return(list(n.cens=n.cens,n.uncens=n.uncens))
}


likelihood <- function(x){
  x$S[x$S==0] <- 1
  x$pi[x$pi==0] <- 1
  likelihood <- sum(x$n.uncens*log(x$pi)+x$n.cens*log(x$S))
  return(likelihood)
}

brier<-function(z,S,time,states,KM,npred){
  (sum(S[,z]^2*(time<=z & states==2)/KM[time]+(1-S[z])^2*(time>z)/KM[z],na.rm=TRUE))/npred
}

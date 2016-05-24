frac.ranks <-
function(x,w=NULL){
  if(is.null(w)) w <- rep(1,length(x))
  x <- sort(x)
  x.unique <- unique(x)
  w.unique <- rep(NA,length(x.unique))
  pik <- rep(NA,length(x.unique))
  Fk  <- rep(NA,length(x.unique))
  k<- 1
  pik[k]    <- weighted.mean(x==x.unique[k],w=w)
  Fk[k] <- 0.5*pik[k]
  for(k in 2:length(x.unique)){
    pik[k]    <- weighted.mean(x==x.unique[k],w=w)
    Fk[k]     <- t(pik[1:(k-1)])%*%rep(1,length(1:(k-1)))+0.5*pik[k]
  }
  Fi <- rep(NA,length(x))
  k <- 1
  for(i in 1:length(x)){
    if(x[i]!=x.unique[k]) k <- k+1
    Fi[i] <- Fk[k]
  }
  return(Fi)
}

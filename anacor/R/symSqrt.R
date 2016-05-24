`symSqrt` <-
function(cc,inv=FALSE)
{
  e<-eigen(cc)
  ev<-e$values
  ind<-which(ev>0)
  lbd<-rep(0,length(ev))
  if (inv) lbd[ind]<-1/sqrt(ev[ind]) else lbd[ind]<-sqrt(ev[ind])
  kk<-e$vectors
  return(kk%*%(lbd*t(kk)))
}


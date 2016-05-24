`fnorm2` <-
function(dmax2,data,alterInt){

  if(any(is.na(dmax2))) return(NaN)


  N <- sum(data)
  m <- colSums(data)
  
  p<-m/N
  p0<-p[1]
  p2<-p[3]
  ro.0.1<-p0*p2/{sqrt(p0*(1-p0))*sqrt(p2*(1-p2))}
  cov<-matrix(c(1,ro.0.1,ro.0.1,1),2,2,byrow=T)
  d<-dmax2

  lower <- switch(alterInt+1, rep(-d,2), rep(-Inf,2), rep(-d,2))
  upper <- switch(alterInt+1, rep( d,2), rep(d, 2), rep(Inf,2))

  
  p.norm<-1-as.numeric(sadmvn(lower=lower, upper=upper,mean=rep(0,2),varcov=cov))
  return(p.norm)
}


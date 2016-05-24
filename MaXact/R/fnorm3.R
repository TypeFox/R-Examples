`fnorm3` <-
  function(dmax3, data, alterInt){
  if(any(is.na(dmax3))) return(NaN)

  N <- sum(data)
  m <- colSums(data)


  p<-m/N
  p0<-p[1]
  p1<-p[2]
  p2<-p[3]
  ro.0.0.5<-p0*(p1+2*p2)/{sqrt(p0*(1-p0))*sqrt((p1+2*p2)*p0+(p1+2*p0)*p2)}
  ro.0.1<-p0*p2/{sqrt(p0*(1-p0))*sqrt(p2*(1-p2))}
  ro.0.5.1<-p2*(p1+2*p0)/{sqrt(p2*(1-p2))*sqrt((p1+2*p2)*p0+(p1+2*p0)*p2)}
  cov<-matrix(c(1,ro.0.0.5,ro.0.1,ro.0.0.5,1,ro.0.5.1,ro.0.1,ro.0.5.1,1),3,3,byrow=T)
  d<-dmax3

  lower <- switch(alterInt+1, rep(-d,3), rep(-Inf,3), rep(d,3))
  upper <- switch(alterInt+1, rep( d,3), rep(d, 3), rep(Inf,3))

  p.norm<-1-as.numeric(sadmvn(lower=lower ,upper=upper, mean=rep(0,3),varcov=cov))
  return(p.norm)
}


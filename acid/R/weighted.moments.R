weighted.moments <-
function(x,w8=NULL){
  if(is.null(w8)) w8<-rep(1,length(x))
  w  <- w8/sum(w8)
  n  <- length(x)
  Ex <- t(x)%*%w
  wm <- Ex
  Ex2<- t(x*x)%*%w 
  wsd<- sqrt(1/(1-sum(w^2))*(Ex2-wm^2)) # using Cramer and cov.wt (using unbiased method)
  wtd.sd<-sqrt(wtd.var(x,w8)) #using Hmisc
  Ex3<- t(x^3)%*%w
  wsk<- n^(5/2)/((n-1)*(n-2))*t(w^(3/2))%*%(((x-wm)/wsd)^3)
  m2.stata<- t(w)%*%((x-wm)^2) # see summarize in stata manual p.1825-6
  m3.stata<- t(w)%*%((x-wm)^3) 
  wsk.stata<-m3.stata / m2.stata^(3/2) 
  list(fm=Ex,weighted.mean=wm,sm=Ex2,weighted.sd=wsd,wtd.sd=wtd.sd,tm=Ex3,w.skew.SAS=wsk,w.skew.Stata=wsk.stata) # using formulae from SAS and from Stata
}

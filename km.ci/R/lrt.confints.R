"lrt.confints" <-
function(time,status,t0,alpha=0.05) {
 ftimes<-sort(unique(time[status==1]))
 k<-length(ftimes)
 dvec<-rep(1,k)
 nvec<-dvec
 jmax<-0
 for(i in 1:k) {
  dvec[i]<-sum(time==ftimes[i] & status==1)
  nvec[i]<-sum(time>=ftimes[i])
  if(t0 >= ftimes[i])
   jmax<-i
 }
 nsub<-nvec[1:jmax]
 dsub<-dvec[1:jmax]
 theta<-qchisq(1-alpha,1)

 #handle no deaths specially
 if( jmax == 0 || sum(dsub) == 0 ) 
  return(list(lower=exp(-theta/(2*max(nsub))),upper=1))
 #else 
 lrt<-function(lambda,opar) { 
  nsub<-opar$nsub
  dsub<-opar$dsub
  opar$theta -2*sum(nsub*log(1+lambda/nsub)-(nsub-dsub)*log(1
   + lambda/(nsub-dsub)))
 }
 plambda<-function(lambda,opar) {
  nsub<-opar$nsub
  dsub<-opar$dsub
  prod(1-dsub/(nsub+lambda))
 }
 t1<-sum(1/(nsub-dsub)-1/nsub)
 l1<- -sqrt(theta/t1)
 #the likelihood is undefined if l1 < lbd
 lbd<- dvec[jmax]-nvec[jmax]
 tol2 <- 1/100
 #if( nvec[jmax]-dvec[jmax]+l1 < 0 ) {
 # l1<-(dvec[jmax]-nvec[jmax]-0.001)/10
 #}
 parlist<-list(nsub=nsub,dsub=dsub,theta=theta)
 for(i in 1:10) {
  if( l1 < lbd ) {
   l1<-lbd + tol2
   tol2<- tol2/2
  }
  v1 <- lrt(l1,parlist)
  if(v1 < 0 ) break
  else {
   slope<-(theta - v1)/l1
   l1<-(theta + 0.5 )/slope
  }
 }
 lower<-bisect(lrt,parlist,l1,0)
 l1<- -l1
 for(i in 1:10) {
  v1<-lrt(l1,parlist)
  if(v1 < 0 ) break
  else {
   slope<-(theta - v1)/l1
   l1<-(theta + 0.5 )/slope
  }
 }
 upper<-bisect(lrt,parlist,0,l1)
 lp<-plambda(lower,parlist)
 up<-plambda(upper,parlist)
 return(list(lower=lp,upper=up))
}


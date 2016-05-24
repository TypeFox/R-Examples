Balanced.Initialization<-function(K,y,nplus,nminus=nplus,eps=1e-12){
  n<-2*nplus
  f<-K%*%y
  Iplus<-seq(y)[y>0]
  Iminus<-seq(y)[y<0]
  fmax<-max(f[Iplus])
  fmin<-min(f[Iminus])
  iplus<-Iplus[match(fmax,f[Iplus],0)]
  iminus<-Iminus[match(fmin,f[Iminus],0)]
###this seems to take into account ties
  lambda<- (fmax-fmin)/2
  beta0<-1-fmax/lambda
  alpha0<-beta0*lambda
### package parameters for the left of the start
  alpha00<-c(slope=beta0,intercept=0)
  
  list(Elbow=c(iplus,iminus),lambda=lambda,alpha0=alpha0,alpha00=alpha00,alpha=rep(1,n))
}
  

logspec2cov<-function(ebeta,vbeta,SARIMA=list(sigma2=1),lags,subdivisions=100){
  nb<-length(ebeta)
  acv<-c()
  pb<-txtProgressBar(min=0,max=lags,style= 3)
  for(k in 0:lags){
  setTxtProgressBar(pb,k)

  dumfun<-function(x){
  #g0<-SARIMAspec(SARIMA,x,log=TRUE)
  g0<-log(SARIMAspec(SARIMA,freq=x)$spec)
  B<-basis(x,nb)
  mu<-B%*%ebeta+g0
  s2<-rowSums((B%*%vbeta)*B)
  exp(mu+0.5*s2)*cos(2*pi*k*x)
  #exp(mu)*cos(2*pi*k*x)
  }

  acv[k+1]<-2*pi*integrate(dumfun,0,0.5,rel.tol=1e-2,subdivisions=subdivisions)$value
  #acv[k+1]<-integrate(dumfun,0,0.5,rel.tol=.Machine$double.eps^0.4,subdivisions=900)$value
  }
  close(pb)
  return(acv)
}


#
# logspec2cov<-function(beta,lags){
#   nb<-length(beta)
#   acv<-c()
#   pb<-txtProgressBar(min=0,max=lags,style= 3)
#   for(k in 0:lags){
#   setTxtProgressBar(pb,k)
#   dumfun<-function(x){exp(basis(x,nb)%*%beta)*cos(2*pi*k*x)}
#   acv[k+1]<-2*pi*integrate(dumfun,0,0.5,subdivisions=500)$value
#   #acv[k+1]<-2*quadcc(dumfun,0,0.5,tol=.Machine$double.eps^0.5)
#   }
#   close(pb)
#   return(acv)
# }

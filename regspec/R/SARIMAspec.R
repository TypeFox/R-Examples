SARIMAspec<-function(model,freq,log=FALSE){
  x<-complex(length(freq),modulus=1,argument=-2*pi*freq)
  p<-length(model$ar);q<-length(model$ma)
  sp<-length(model$seasonal$sar);sq<-length(model$seasonal$sma)
  xmat<-outer(x,0:max(p,q),FUN="^")
  
  if(is.null(model$seasonal$period)){model$seasonal$period<-1}
  sxmat<-outer(x,model$seasonal$period*0:max(0,sp,sq),FUN="^")
  
  if(!is.null(model$ar)){model$ar<--model$ar}
  if(!is.null(model$seasonal$sar)){model$seasonal$sar<--model$seasonal$sar}
  
  logspec<-log(model$sigma2)-
    2*log(Mod(xmat[,1+0:p,drop=FALSE]%*%c(1,model$ar)))+
    2*log(Mod(xmat[,1+0:q,drop=FALSE]%*%c(1,model$ma)))-
    2*log(Mod(sxmat[,1+0:sp,drop=FALSE]%*%c(1,model$seasonal$sar)))+
    2*log(Mod(sxmat[,1+0:sq,drop=FALSE]%*%c(1,model$seasonal$sma)))
    if(!is.null(model$d)){logspec<-logspec-2*model$d*log(Mod(1-x))}
    if(!is.null(model$seasonal$sd)){logspec<-logspec-2*model$seasonal$sd*log(Mod(1-x^model$seasonal$period))}
  
  out<-list(freq=freq,logspec=logspec)
  if(log==FALSE){out$logspec<-NULL;out$spec<-exp(logspec)}
  return(out)
}

# wplot<-seq(0,0.5,length=120)
# SARIMA0<-list(ar=c(-0.5, 0.4, 0.8),ma=0.2,sigma2=1,seasonal=list(sar=0.5,period=12))
# #SARIMA0<-list(ar=c(-0.5, 0.4, 0.8),ma=0.2,sigma2=1)
# spec.true<-SARIMAspec(SARIMA0,freq=wplot)
# plot(wplot,log(spec.true$spec),type="l")
# 
# spec.true<-ARMAspec(model=SARIMA0,freq=wplot,plot=FALSE)
# points(wplot,log(spec.true$spec),type="l",col=2)

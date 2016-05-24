####Gaussian Quadrature Part Integration
Innerfunc=function(zz,rstat,model,
                   ww2=0.8049141,ww1=0.08131284,ww3=0.8049141,ww4=0.08131284,
                   xx1=-1.65068,xx2=-0.5246476,xx3= 0.5246476,xx4= 1.65068){
  if(zz==1){
    ff1=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx1-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    ff2=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx2-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    ff3=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx3-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    ff4=(1+0.5*((rstat-sqrt(2)*model$parameter$eta1*xx4-model$parameter$gamma1)^2)/model$parameter$beta1)^(-model$parameter$alpha1-0.5)
    
    outvalue=exp(lgamma(model$parameter$alpha1+0.5)-lgamma(model$parameter$alpha1))/sqrt(pi*model$parameter$beta1)
    
  }else{
    ff1=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx1-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    ff2=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx2-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    ff3=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx3-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    ff4=(1+0.5*((rstat-sqrt(2)*model$parameter$eta0*xx4-model$parameter$gamma0)^2)/model$parameter$beta0)^(-model$parameter$alpha0-0.5)
    
    outvalue=exp(lgamma(model$parameter$alpha0+0.5)-lgamma(model$parameter$alpha0))/sqrt(pi*model$parameter$beta0)
    
  }
  ff<-c(ff1,ff2,ff3,ff4)
  integralvalue=outvalue*(c(ww1,ww2,ww3,ww4)%*%ff)
  return(integralvalue)
}

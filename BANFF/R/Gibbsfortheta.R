#####giibs for mu ars for var
Gibbsfortheta<-function(initalmu,initalvar,index,model)
{
  uprr=initalmu
  initialvarvar=initalvar
  
  if(index==0){
    
    
    myvar<-rigamma(1,model$parameter$alpha0,model$parameter$beta0)
    
    
    meanforupmu=(myvar*model$parameter$gamma0+model$parameter$eta0^2*uprr)/(myvar+model$parameter$eta0^2)
    varforupmu=myvar*model$parameter$eta0^2/(myvar+model$parameter$eta0^2)
    
    
    mymu<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
    
    
  }else{
    
    
    myvar<-rigamma(1,model$parameter$alpha1,model$parameter$beta1)
    
    
    
    meanforupmu=(myvar*model$parameter$gamma1+model$parameter$eta1^2*uprr)/(myvar+model$parameter$eta1^2)
    varforupmu=myvar*model$parameter$eta1^2/(myvar+model$parameter$eta1^2)
    
    mymu<-rnorm(1,mean=meanforupmu,sd=sqrt(varforupmu))
    
    
  }
  
  
  return(c(mymu,myvar))
}


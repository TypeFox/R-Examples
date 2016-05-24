# Density
charact.MTSgam <- function(t,mu0,mu,sig,a,alpha,lambda_p,lambda_m,
                           Parametrization){
  if(alpha==2){
    if(Parametrization=="A"){
      ycgf<-1i*t*mu0-a*log(1-sig^2*(mu*1i*t-t^2/2))
    } 
    if(Parametrization=="B"){
      ycgf<-1i*t*mu0-a*log(1-(mu*1i*t-sig^2*t^2/2))
    }
  }else{
    L_cts<-(((lambda_p-1i*t)^alpha-lambda_p^alpha+(lambda_m+1i*t)^alpha-lambda_m^alpha)/(alpha*(alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2)))
          +(1i*t*(lambda_p^(alpha-1)-lambda_m^(alpha-1)))/((alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2))))
    if(Parametrization=="A"){
    #  L_cts_y<-sig^2*L_cts
      ycgf<-1i*t*mu0-a*log(1-sig^2*(mu*1i*t+L_cts))
    }
    if(Parametrization=="B"){
      L_cts_y<-sig^2*L_cts
      ycgf<-1i*t*mu0-a*log(1-(mu*1i*t+L_cts_y))
    }
  }
  y<-exp(ycgf)
  return(y)
}

charact.User<-function(t,mu0,mu,sig,
             a,alpha,lambda_p,lambda_m,UseMGF,paramMixing,Parametrization){
  if(alpha==2){
    if(Parametrization=="A"){
      ySub<-sig^2*(mu*1i*t+t^2/2)
    }
    if(Parametrization=="B"){
      ySub<-mu*1i*t+sig^2*t^2/2
    }
  }else{
    L_cts<-(((lambda_p-1i*t)^alpha-lambda_p^alpha+(lambda_m+1i*t)^alpha-lambda_m^alpha)/(alpha*(alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2)))
    +(1i*t*(lambda_p^(alpha-1)-lambda_m^(alpha-1)))/((alpha-1)*(lambda_p^(alpha-2)+lambda_m^(alpha-2))))
    if(Parametrization=="A"){
      ySub<-sig^2*(mu*1i*t + L_cts)
    }
    if(Parametrization=="B"){
      ySub<-mu*1i*t+sig^2*L_cts    #We apply the result in Edit Mercuri 2014. The Charatacteristic of MixedTS can be written as time-changed Levy process
    }
  }
  ycgf<-1i*t*mu0+UseMGF(u=ySub,paramMixing=paramMixing,a=a)
  y<-exp(ycgf)
  return(y)
}


ChFunToDens.MTSgam <- function(n, alim, blim, mu0, mu, sig, a, alpha,
                               lambda_p,lambda_m, Parametrization) {
  # Internal function for transforming the characteristic function in Density,
  # using fft
  i <- 0:(n-1)            
  dx <- (blim-alim)/n     
  x <- alim + i * dx      
  dt <- 2*pi / ( n * dx ) 
  c <- -n/2 * dt          
  d <-  n/2 * dt          
  t <- c + i * dt         
  phi_t <- charact.MTSgam(t,mu0,mu,sig,
                           a,alpha,lambda_p,
                          lambda_m, Parametrization)
  X <- exp( -(0+1i) * i * dt * alim ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
}

ChFunToDens.User<-function(n,alim,blim,
                 mu0,mu,sig,a,alpha,lambda_p,
                 lambda_m,UseMGF,paramMixing,Parametrization){

  i <- 0:(n-1)            
  dx <- (blim-alim)/n     
  x <- alim + i * dx      
  dt <- 2*pi / ( n * dx ) 
  c <- -n/2 * dt          
  d <-  n/2 * dt          
  t <- c + i * dt         
  phi_t <- charact.User(t,mu0,mu,sig,
                          a,alpha,lambda_p,lambda_m,UseMGF,paramMixing,
                        Parametrization)
  X <- exp( -(0+1i) * i * dt * alim ) * phi_t
  Y <- fft(X)
  density <- dt / (2*pi) * exp( - (0+1i) * c * x ) * Y
  data.frame(
    i = i,
    t = t,
    characteristic_function = phi_t,
    x = x,
    density = Re(density)
  )
  
}

dMixedTS.aux <- function(x,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,MixingDens="Gamma",N=2^10,
                         UseMGF=NULL, paramMixing=NULL, Parametrization = "A"){
# Density function of a 
# Mixed Tempered Stable 
# distribution developed in [Rroji and Mercuri 2014]
  if(length(x)==1){
    alim<-min((-abs(x)-0.5),setInf)
    blim<-max((abs(x)+0.5),setSup)
  }else{
    xdummy<-na.omit(x[is.finite(x)])
    alim<-min(min(xdummy)-0.1,setInf)
    blim<-max(max(xdummy)+0.1,setSup)
  }
  if(MixingDens=="Gamma"){
    
    invFFT<-ChFunToDens.MTSgam(n=N,alim=alim,blim=blim,mu0=mu0,mu=mu,sig=sig,a=a,alpha=alpha,lambda_p=lambda_p,lambda_m=lambda_m, 
                               Parametrization=Parametrization)
    dens<-approx(invFFT$x,invFFT$density,x)
    return(dens$y)
  }
  if(MixingDens=="User"){
    invFFT<-ChFunToDens.User(n=N,alim=alim,blim=blim,
                             mu0=mu0,mu=mu,sig=sig,a=a,
                             alpha=alpha,lambda_p=lambda_p,
                             lambda_m=lambda_m,UseMGF=UseMGF,
                             paramMixing=paramMixing, 
                             Parametrization=Parametrization)
    dens<-approx(invFFT$x,invFFT$density,x)
    return(dens$y)
  }
}

dMixedTS.saddlepoint <- function(x,mu0,mu,sig,a,alpha,
                                 lambda_p,lambda_m,setSup=NULL,
                                 setInf=NULL,MixingDens="Gamma", 
                                 Parametrization=Parametrization){
    
  
}

pMixedTS.aux <- function(q,mu0,mu,sig,a,alpha,lambda_p,lambda_m,lower=-Inf,
                     setSup=NULL,setInf=NULL,MixingDens="Gamma",
                     UseMGF=NULL, paramMixing=NULL,N, 
                     Parametrization="A"){
  
  env <- new.env()
  env$mu0 <- mu0
  env$mu <- mu
  env$a <- a
  env$sig <- sig
  env$alpha <- alpha
  env$lambda_p <- lambda_p
  env$lambda_m <- lambda_m
  env$setSup <- setSup
  env$setInf <- setInf
  env$MixingDens <- MixingDens
  env$UseMGF<-UseMGF
  env$paramMixing<-paramMixing
  env$N<-N
  env$Parametrization<-Parametrization
 # if(MixingDens=="Gamma"){ 
          myfun<- function(x,env){
            env$mu0 -> mu0
            env$mu -> mu
            env$a -> a
            env$sig -> sig
            env$alpha -> alpha
            env$lambda_p -> lambda_p
            env$lambda_m -> lambda_m
            env$setSup -> setSup
            env$setInf -> setInf
            env$MixingDens -> MixingDens
            env$UseMGF->UseMGF
            env$paramMixing->paramMixing
            env$N->N1
            env$Parametrization -> Parametrization
            n<-length(x)
            y<-dMixedTS.aux(x,mu0=mu0,mu=mu,sig=sig,a=a,
                            alpha=alpha,lambda_p=lambda_p,
                            lambda_m=lambda_m,setSup=setSup,
                            setInf=setInf,MixingDens=MixingDens,
                            N=N1,UseMGF=UseMGF, 
                            paramMixing=paramMixing,
                            Parametrization=Parametrization)
   }
    lower<--Inf
      if(!is.finite(lower)){
        #lower<-(min(q[is.finite(q)])-1)
        lower<--20
        if(!is.null(setInf)){
          lower<-setInf-5
        }
      }
    approach<-"RectangularMethod"
    if(approach=="RectangularMethod"){
      qmax<-max(q[is.finite(q)])
      NPartition<-12
      point<-seq(lower,qmax,length=2^NPartition)
      lengthStep<-diff(point)[1]
      midPoint<-lengthStep/2+point
      probpoint<-cumsum(myfun(x=midPoint[-length(midPoint)],env=env))*lengthStep
      prob<-approx(point[-1],probpoint,q)
      # we compute the integral between -\infty to lower
#       Udummy<-seq(0,1,length=2^NPartition)
#       
#       ldummy<-diff(Udummy)[1]
#       stepU<-(Udummy+ldummy/2)
#       midDummy<-lower-(1-(stepU))/(stepU)
#       fin<-length(midDummy)
#       y<-myfun(x=midDummy[-fin],env=env)/stepU[-fin]^2
#       y[is.na(y)]=0
#       probdummy<-sum(y)*ldummy
#       prob$y<-as.numeric(prob$y)+probdummy
       return(prob$y)
    }else{
      if(approach=="Simpson"){
        warning("We will implement as son as possible!")
        return(NULL)
      }  
    }
  
}



# Quantile function

qMixedTS.aux <- function(p,mu0,mu,sig,a,
                         alpha,lambda_p,lambda_m,
                         setSup=NULL,setInf=NULL,MixingDens="Gamma",
                         N=N,UseMGF, paramMixing,Parametrization="A",
                         interval=c(-100,100)){
  # It is enough to invert the function pMixedTS.aux using the function \texttt{uniroot} 
  # of the stats package
  env <- new.env()
  env$mu0 <- mu0
  env$mu <- mu
  env$a <- a
  env$sig <- sig
  env$alpha <- alpha
  env$lambda_p <- lambda_p
  env$lambda_m <- lambda_m
  env$setSup <- setSup
  env$setInf <- setInf
  env$MixingDens <- MixingDens
  env$UseMGF<-UseMGF
  env$paramMixing<-paramMixing
  env$N<-N
  env$Parametrization<-Parametrization
  myObjectFun<-function(q,prob,env){

    env$mu0 -> mu0
    env$mu -> mu
    env$a -> a
    env$sig -> sig
    env$alpha -> alpha
    env$lambda_p -> lambda_p
    env$lambda_m -> lambda_m
    env$setSup -> setSup
    env$setInf -> setInf
    env$MixingDens -> MixingDens
    env$UseMGF->UseMGF
    env$paramMixing->paramMixing
    env$N->N1
    env$Parametrization->Parametrization
    y <- pMixedTS.aux(q=q,mu0=mu0,mu=mu,sig=sig,a=a,alpha=alpha,
                  lambda_p=lambda_p,lambda_m=lambda_m,setSup=setSup,
                  setInf=setInf,MixingDens=MixingDens,
                  UseMGF=UseMGF,paramMixing=paramMixing,
                  Parametrization=Parametrization,
                  N=N1)-prob
  }
  if(length(p)==1){
    res <-uniroot(f=myObjectFun,interval=interval,prob=p,env=env)$root
#    res <- fsolve(f=myObjectFun,x0=0.1,prob=p,env=env)$x
    return(res)
  }else{
    probmax <- max(p)
    probmin <- min(p)
    if(!is.null(setSup)&&!is.null(setInf)){
      interval<-c(setInf,setSup)
    }
    maxres <- uniroot(f=myObjectFun,interval=interval,prob=probmax,env=env)$root
    minres <- uniroot(f=myObjectFun,interval=interval,prob=probmin,env=env)$root
    stepquant<-seq(minres,maxres,length=2^12)
    probabilityTS<-pMixedTS.aux(stepquant,mu0,mu,sig,a,alpha,lambda_p,lambda_m,lower=-Inf,
                         setSup=setSup,setInf=setInf,MixingDens=MixingDens,
                         UseMGF=UseMGF,paramMixing=paramMixing,N=env$N,Parametrization=Parametrization)
    res<-approx(probabilityTS,stepquant,p)
    return(res$y)
  }
  
}

# Random Number generetor

rMixedTS.aux <- function(n,mu0,mu,sig,a,alpha,lambda_p,lambda_m,setSup=NULL,setInf=NULL,
                         MixingDens="Gamma",UseMGF, paramMixing, Parametrization="A", Nstep=2^12){
  UniformNumber  <- runif(n, min = 0.0001, max = 0.999)
  MixedTSNumber  <- qMixedTS.aux(UniformNumber,mu0,mu,sig,a,
                                 alpha,lambda_p,lambda_m,
                                 setSup=setSup,setInf=setInf,
                                 MixingDens=MixingDens, 
                                 Parametrization=Parametrization, 
                                 interval=c(-100,100),
                                 UseMGF=UseMGF, 
                                 paramMixing=paramMixing, 
                                 N=Nstep)
return(MixedTSNumber)
}

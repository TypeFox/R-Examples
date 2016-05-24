WLP.test<-function(x,type="CQ",funEstimate="MLE", s1=-1,s2=0.4,mr=100,a=-5,nsim=500){

  ##Family of the test statistics based on the Laplace transform
  LP.statistic<-function(x,type,funEstimate,s1,s2,mr,a){
    TYPE <- deparse(substitute(type))
    funEst <- deparse(substitute(funEstimate))
    if(funEstimate=="MLE") {y.est<-MLEst(x)}
    else if (funEstimate=="ME"){y.est<-MEst(x)}
    else if (funEstimate=="LSE"){y.est<-LSEst(x)}
    else{stop(paste("unknown estimation method ", funEst, "!"))}
    y<-y.est$y
    z <- exp(-exp(y)) 
    n<-length(x)
    l <- seq(length=n,from=1,to=n)
    d_2gamma <-function(x){   
      d_2gamma = (trigamma(x)+digamma(x)*digamma(x))*gamma(x) 
    }
    
    ###### Fonction that computes the Cov between to values (ME) 
    CovME <- function(s1,s2){
      Cov=0
      #constante d'euler
      eul=-digamma(1)
      #skewness ev(0,1)
      a3=(12*sqrt(6)/pi^3)*1.202057  #zeta(3)=1.202057
      #kurtosis
      a4=5.4
      d1_s1=digamma(1-s1)*gamma(1-s1)
      d1_s2=digamma(1-s2)*gamma(1-s2)
      d2_s1=d_2gamma(1-s1)
      d2_s2=d_2gamma(1-s2)
      Cov = gamma(1-s1-s2)-gamma(1-s1)*gamma(1-s2)
      Cov = Cov +s2*gamma(1-s2)*(d1_s1+eul*gamma(1-s1))+3*s2/pi^2*(eul*gamma(1-s2)+d1_s2)*(d2_s1+2*eul*d1_s1+(eul^2-pi^2/6)*gamma(1-s1))
      Cov = Cov + s1*gamma(1-s1)*(d1_s2+eul*gamma(1-s2))+3*s1/pi^2*(eul*gamma(1-s1)+d1_s1)*(d2_s2+2*eul*d1_s2+(eul^2-pi^2/6)*gamma(1-s2))
      Cov = Cov -s2*s1*gamma(1-s2)*(-pi^2/6*gamma(1-s1)+eul*a3*pi/(2*sqrt(6))*gamma(1-s1)+pi*a3/(2*sqrt(6))*d1_s1)   
      CovME <- Cov +3*s1*s2/pi^2*(eul*gamma(1-s2)+d1_s2)*(-a3*pi^3/(6*sqrt(6))*gamma(1-s1)+pi^2/12*(a4-1)*(eul*gamma(1-s1)+d1_s1))   
    }
    
    ###### Fonction that computes the Cov between to values (MLE) 
    CovMLE <- function(s1=-0.1,s2=0.02){
      Cov=0
      #constante d'euler
      eul=-digamma(1)
      d1_s1=digamma(1-s1)*gamma(1-s1)
      d1_s2=digamma(1-s2)*gamma(1-s2)
      d2_s1=digamma(2-s1)*gamma(2-s1)
      d2_s2=digamma(2-s2)*gamma(2-s2)
      Cov=gamma(1-s1-s2)-gamma(1-s1)*gamma(1-s2)
      ds1=NULL
      ds2=NULL
      ds1=c(-s1*gamma(1-s1))
      ds1=c(ds1,s1*d1_s1)   
      ds2=c(-s2*gamma(1-s2))
      ds2=c(ds2,s2*d1_s2)
      a=1+6/pi^2*(1-eul)^2
      b=-6/pi^2*(1-eul)
      A = matrix(data=c(a,b,b,6/pi^2),nrow=2,ncol=2)
      b=eul-1
      l1=NULL
      l2=NULL
      a1=-gamma(1-s1)-d1_s1+d2_s1
      a2=-gamma(1-s2)-d1_s2+d2_s2
      l1=c(-gamma(1-s1)+gamma(2-s1),a1)
      l2=c(-gamma(1-s2)+gamma(2-s2),a2)
      Cov=Cov+ds2%*%A%*%l1
      Cov=Cov+ds1%*%A%*%l2
      CovMLE <- Cov+ds1%*%t(A)%*%ds2 
    }
    
    #### Test statistic of Cabana and Quiroz
    if(type=="CQ"){
      if(funEstimate=="MLE"){Cov=CovMLE}
      else if(funEstimate=="ME"){Cov=CovME}
      else {stop(paste("statistic CQ is not defined in the case of","",funEst,"!"))}
      e1 = exp(-y*s1)
      e2 = exp(-y*s2)
      v1 = Cov(s1,s1)
      v2 = Cov(s2,s2)
      v12 = Cov(s1,s2)
      vn1 = sqrt(n)*(mean(e1)-gamma(1-s1))
      vn2 = sqrt(n)*(mean(e2)-gamma(1-s2))
      CQ = 1/(v1*v2-v12^2)
      LP.statistic <- CQ*(v2*vn1^2-2*vn1*vn2*v12+v1*vn2^2)
      ####Test statistic based on the Laplace transform LT
    }else if(type=="LT"){
      l=seq(-mr,-1)
      l=l/mr
      m=matrix(data=l,nrow=length(l),ncol=n)
      e=exp(-y*t(m))
      s=colSums(e)

      LP.statistic <- n*sum((exp(-exp(a*l)+a*l)*(gamma(1-l)-1/n*s)^2))  
     ########## CQ^*
    }else if(type=="CQ*"){
     s1 <- -0.1
     s2 <- 0.02
     e1 = exp(-y*s1)
     e2 = exp(-y*s2)
     v1 = 1.59
     v2 = 0.53
     v12 = 0.91
     vn1 = sqrt(n)*(mean(e1)-gamma(1-s1))
     vn2 = sqrt(n)*(mean(e2)-gamma(1-s2))
     Qn = 1/(v1*v2-v12^2)
     LP.statistic <- Qn*(v2*(vn1)^2-2*vn1*vn2*v12+v1*(vn2)^2)}

    return(list(statistic=as.vector(LP.statistic),eta=y.est$eta,beta=y.est$beta))
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  EST <- NULL
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="CQ"){
    METHOD="Smooth test statistic of Cabana and Quiroz for the Weibull distribution"
  }else if(as.character(type)=="LT"){
    METHOD="Smooth test statistic based on the Laplace transform for the Weibull distribution"
  }else if(as.character(type)=="CQ*"){
    METHOD="Smooth test statistic CQ* for the Weibull distribution"}
  else  stop(paste("The chosen method ",TYPE," is unknown"))
  
  EST <- deparse(substitute(funEstimate))
  if(as.character(funEstimate)=="MLE"){EST="using the MLEs "
  } else if(as.character(funEstimate)=="ME"){EST="using the MEs "
  } else if(as.character(funEstimate)=="LSE"){EST="using the LSEs "
  } else  stop(paste("The chosen estimation method ",EST," is unknown"))
  stat <- LP.statistic(x,type,funEstimate,s1,s2,mr,a)
  statistic.obs=stat$statistic
  estimate.obs <- c(stat$eta,stat$beta)
  fun<-function(y){
    fun <- LP.statistic(y,type,funEstimate,s1,s2,mr,a)
    return(fun$statistic)
  }
  sim.statistic <- GoFsim(nsim,n,fun)
  p_val <- sum(sim.statistic>=statistic.obs)/nsim

  WLP.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                   method=paste(METHOD, " ", EST),
                   estimate = c(eta=estimate.obs[1],beta=estimate.obs[2]),  
                   data.name=DNAME)
  
  class(WLP.test) <-"htest"
  return(WLP.test)
}

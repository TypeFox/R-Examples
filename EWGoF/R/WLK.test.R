WLK.test<-function(x,type="GG1",funEstimate="MLE",procedure="S",nsim=500,r=0){
  ##Family of the Likelihood based test
  LK.statistic<-function(x,type,funEstimate,procedure,r=0){
    TYPE <- deparse(substitute(type))
    PROC <- deparse(substitute(procedure))
    funEst <- deparse(substitute(funEstimate))
    if(funEstimate=="MLE") {
      if(r==0){y.est<-MLEst(x)}
      else y.est<-MLEst_c(x,r)
    }
    else if (funEstimate=="ME"){y.est<-MEst(x)}
    else if (funEstimate=="LSE"){y.est<-LSEst(x)}
    #else if (funEstimate=="BLOM"){y.est<-BLOMEst(x)}
    else{stop(paste("unknown estimation method ", funEst, "!"))} 
    y <- sort(y.est$y)
    n <- length(x)
    #### Test statistic based on the Generalized Gamma
    if(type=="GG1"){
      ##Function that approximates the inverse of the digamma function
      sol_minka<-function(k){
        if(k< -2.22){ sol_minka = -1/(k-digamma(1))
        }else{ sol_minka = exp(k)+.5}
        for(i in 1:4){
          t=sol_minka
          sol_minka=sol_minka-(digamma(sol_minka)-k)/trigamma(sol_minka)}
        s=sol_minka} 
      #y.est <- MLEst(x)
      y <- sort(y.est$y)
      k <- sol_minka(mean(y))
      A <- sum(y)-n*digamma(1)
      B <- n*trigamma(1)
      S <- A^2/B
      W <-  B*(1-k)^2
      LR <- 2*(-n*log(gamma(k))+(k-1)*sum(y))
      
      #### Test statistic based on the Generalized Gamma after transformation GG2
    }else if(type=="GG2"){
      #Function to estimate k hat (the MLE) of the Generalized Gamma distribution
      f1 <- function(toto,vect){        
        x=exp(vect/sqrt(toto))
        f1=n*toto*log(toto)+n*(toto-.5)-n*toto*digamma(toto)
        f1=f1+.5*sqrt(toto)*sum(vect)-toto*sum(x)+0.5*sqrt(toto)*sum(vect*x)
        f1=(f1)^2
      }     
      lnt <- optimize(f1,c(0.0001,4),maximum=FALSE,vect=y,tol=1/10^5)
      k <- lnt$minimum
      l <- exp(y/sqrt(k))
      A <- n/2-n*digamma(1)+.5*sum(y)-sum(exp(y))+.5*sum(y*exp(y))
      B <- -(3*n/2-n*trigamma(1)-.25*sum(y)+0.25*sum(y*exp(y))-.25*sum((y^2)*exp(y)))
      S <- A^2/B
      W <- (k-1)^2*B
      LR <- n*(k-.5)*log(k)-n*log(gamma(k))+(sqrt(k)-1)*sum(y)
      LR <- LR+sum(exp(y))-k*sum(l)
      LR<- 2*LR
      ####Test statistic based on the Exponentiated Weibull
    }else if(type=="EW"){
      #MLE of the third parameter alpha of the distribution 
      alpha <- -n/sum(log(1-exp(-exp(y))))
      W <- n*(1-alpha)^2
      S <- n*(1-1/alpha)^2
      LR <- 2*(n*log(alpha) - n + n/alpha)
      ####Test statistic based on the Power Generalized Weibull   
    }else if(type=="PGW"){
      #MLE of the third parameter lambda
      f1 <- function(toto,vect){
        a=(1+exp(vect))^(1/toto)
        f1=n*toto+sum(log(1+exp(vect)),na.rm = FALSE)
        f1=f1-sum(log(1+exp(vect))*a,na.rm = FALSE)
        f1=abs(f1)
      }
      lnt <- optimize(f1,c(0.0001,10),maximum=FALSE,vect=y,tol=1/10^5)
      lambda <- lnt$minimum
      t <- (1+exp(y))
      A <- n - sum(log(t)*exp(y))
      B <- n - 2*sum(log(t)*(t-1))-sum(t*(log(t))^2)
      S <- -A^2/B
      W <- -B*(1-lambda)^2
      LR <- 2*(-n*log(lambda)+(1/lambda-1)*sum(log(t))-sum(t^(1/lambda))+sum(t))
      
      ### Test statistic based on Marshall-Olkin distribution
    }else if(type=="MO"){
      
      #MLE of the third parameter
      f2 <- function(toto,vect){
        v = 1-(1-toto)*exp(-exp(vect))
        f2 = length(vect)-2*toto*sum(exp(-exp(vect))/v,na.rm = FALSE)
        f2 = abs(f2)
      }
      lnt<-optimize(f2,c(0.000001,10),maximum=FALSE,vect=y,tol=1/10^5)
      mu <- lnt$minimum
      A <- n-2*sum(exp(-exp(y)))
      B <- n-2*sum(exp(-2*exp(y)))
      S <- A^2/B
      W <- B*(mu-1)^2
      LR <- 2*(n*log(mu)-sum(log((1-(1-mu)*exp(-exp(y)))^2)))
      ###Test statistic based on the Modified Weibull distribution
    }else if(type=="MW"){
      #MLE of the third parameter
      f1 <- function(toto,vect){
        z = toto*exp(vect)
        v = 1/(exp(-vect) + toto)
        f1 = sum(exp(vect)) + sum(v) - sum(exp(z + 2*vect))
        f1 = abs(f1)
      }
      lnt <- optimize(f1,c(0.00001,10),maximum=FALSE,vect=y,tol=1/10^5)
      lam <- lnt$minimum
      s <- 1/(exp(-y)+lam)^2
      A <- 2*sum(exp(y))-sum(exp(2*y))
      B <- sum(exp(2*y)+exp(3*y))
      W <- lam^2*B
      S <- A^2/B
      LR <- 2*((lam+1)*sum(exp(y))+sum(log(1+lam*exp(y))))
      LR<- LR-2*sum(exp(y)*exp(lam*exp(y)))  
    }
    if(procedure=="S"){LK.statistic <- S}
    else if(procedure=="W"){LK.statistic <- W}
    else if(procedure=="LR"){LK.statistic <- LR}
    else{stop(paste(PROC, "does not exist choose S, W or LR"))}
    
    return(LK.statistic=list(statistic=LK.statistic,eta=y.est$eta,beta=y.est$beta))
  }
  
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  EST <- NULL
  PROC <- NULL
  family<-1
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="GG1"){
    METHOD="Test based on the Generalized Gamma distribution for the Weibull distribution"
  }else if(as.character(type)=="GG2"){
    METHOD="Test based on the Generalized Gamma distribution (2) for the Weibull distribution "
  }else if(as.character(type)=="EW"){
    METHOD="Test based on the Exponentiated Weibull distribution for the Weibull distribution "
  }else if(as.character(type)=="PGW"){
    METHOD="Test based on the Power Generalized Weibull distribution for the Weibull distribution "
  }else if(as.character(type)=="MO"){
    METHOD="Test based on the Marshall Olkin distribution for the Weibull distribution "
  }else if(as.character(type)=="MW"){
    METHOD="Test based on the Modified Weibull distribution for the Weibull distribution "
  }else if(as.character(type)=="T"){
    METHOD="Test based on the combination of MW and PGW for the Weibull distribution "
    family=2
  }else if(as.character(type)=="G"){
    METHOD="Test based on the combination of GG statistic "
    family=2}
  else  stop(paste("The chosen method ",TYPE," is unknown"))
  if(family==1){
    if(as.character(type)=="GG1"&& as.character(funEstimate)!="MLE"){
    print("GG1 are only defined with the MLEs!")
    funEstimate="MLE"
    EST="using the MLEs"}
    EST <- deparse(substitute(funEstimate))
  if(as.character(funEstimate)=="MLE"){EST="using the MLEs "
  } else if(as.character(funEstimate)=="ME"){EST="using the MEs "
  } else if(as.character(funEstimate)=="LSE"){EST="using the LSEs "
  } else  stop(paste("The chosen estimation method ",EST," is unknown"))
  if(as.character(procedure)=="S"){PROC="with the Score procedure "
  } else if(as.character(procedure)=="W"){PROC="with Wald procedure"
  } else if(as.character(procedure)=="LR"){PROC="with the Likelihood ratio procedure"
  } else  stop(paste("The chosen procedure ",PROC," is unknown"))  
  stat <- LK.statistic(x,type,funEstimate,procedure,r)
  statistic.obs=stat$statistic
  estimate.obs <- c(stat$eta,stat$beta)
  fun<-function(y){
    fun <- LK.statistic(y,type,funEstimate,procedure,r)
    return(fun$statistic)
  }
 }else if(as.character(type)=="T"){
    stat1 <- LK.statistic(x,"MW","MLE","W")
    statistic.obs1=stat1$statistic
    stat2 <- LK.statistic(x,"PGW","ME","W")
    statistic.obs2=stat2$statistic
    estimate.obs <- c(stat1$eta,stat1$beta)
    fun1<-function(y){
     fun1 <- LK.statistic(y,"MW","MLE","W")
     return(fun1$statistic)
    }
    fun2<-function(y){
    fun2 <- LK.statistic(y,"PGW","ME","W")
    return(fun2$statistic)
  }
  sim <- GoFsim2d(nsim,n,fun1,fun2)   
  sim.statistic1 <- sim[1,]
  sim.statistic2 <- sim[2,]
  statistic.obs <- 0.5*abs(statistic.obs1-mean(sim.statistic1))/sd(sim.statistic1)
  statistic.obs <- statistic.obs + 0.5*abs(statistic.obs2-mean(sim.statistic2))/sd(sim.statistic2)
  sim.statistic <- 0.5*abs((sim.statistic1-mean(sim.statistic1)))/sd(sim.statistic1)
  sim.statistic <- sim.statistic + 0.5*abs((sim.statistic2-mean(sim.statistic2)))/sd(sim.statistic2)
 }else{
    stat1 <- LK.statistic(x,"GG1","MLE","W",r)
    statistic.obs1=stat1$statistic
    stat2 <- LK.statistic(x,"GG1","MLE","LR",r)
    statistic.obs2=stat2$statistic
    estimate.obs <- c(stat1$eta,stat1$beta)
    fun1<-function(y){
      fun1 <- LK.statistic(y,"GG1","MLE","W",r)
      return(fun1$statistic)
    }
    fun2<-function(y){
      fun2 <- LK.statistic(y,"GG1","MLE","LR",r)
      return(fun2$statistic)
    }
    sim <- GoFsim2d(nsim,n,fun1,fun2)   
    sim.statistic1 <- sim[1,]
    sim.statistic2 <- sim[2,]
    statistic.obs <- 0.5*abs(statistic.obs1-mean(sim.statistic1))/sd(sim.statistic1)
    statistic.obs <- statistic.obs + 0.5*abs(statistic.obs2-mean(sim.statistic2))/sd(sim.statistic2)
    sim.statistic <- 0.5*abs((sim.statistic1-mean(sim.statistic1)))/sd(sim.statistic1)
    sim.statistic <- sim.statistic + 0.5*abs((sim.statistic2-mean(sim.statistic2)))/sd(sim.statistic2)  
   }

  sim.statistic<-switch(type,"T" = sim.statistic,"G" = sim.statistic, GoFsim(nsim,n,fun))
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  WLK.test <- list(statistic =c(S= statistic.obs), p.value = p_val, 
                 method = paste(METHOD," ", EST, " ", PROC),
                 estimate = c(eta = estimate.obs[1], beta = estimate.obs[2]),  
                 data.name = DNAME)
  
  class(WLK.test) <-"htest"
  return(WLK.test)
}


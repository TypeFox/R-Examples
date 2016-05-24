
LRI.test<-function(x,type="BH", a=1, nsim=2000){
  #Family of the test statistics based on the Laplace, the integrated distribution 
  #function and the mean residual life 
  LRI.statistic<-function(x,type="BH",a=1){
    TYPE <- deparse(substitute(type))
    n = length(x)
    #EMV
    e = n/sum(x)    
    y=sort(e*x)
    ### Test of Baringhaus and Henze
    if(type=="BH"){
      LRI.statistic <- GoFBH(x,a)
    }
    ### Test of Henze
    else if(type=="He"){
      a1 <- seq(1,n)
      integrand <- function(y) {exp(-y)/(y)}
      a2 <- integrate(integrand, lower = a, upper = Inf)$value
      he <- n*(1-a*exp(a)*a2)  
      s <- he
      for (j in 1:n){
        t=y[j]+a
        a1[j]=integrate(integrand, lower = t, upper = Inf)$value
      }
      LRI.statistic <-s+GoFHe(x,a,a1)
  
    } 
    ### Test of Klar
    else if(type=="Kl"){
      k<-1:n
      s<- GoFKl(x,a)
      b <- s[1]
      c <- s[2]
      if (a==0){
        LRI.statistic<-n/2-2*sum(exp(-y))-1/(3*n)*(t(n-k-1)%*%y^3)+n^(-1)*b 
      } else { 
        LRI.statistic<-2*(3*a+2)*n/((2+a)*(1+a)^2)-2*a^3/(1+a)^2*sum(exp(-(1+a)*y))-2/n*sum(exp(-a*y))+2/n*c 
        
      }
    }
    ### Test of Baringhaus based on the integrated distribution function (Cramer-von mises)
    else if(type=="BHC"){
      LRI.statistic <-GoFBHC(x,a)
    }
    ###Test of Baringhaus based on the integrated distribution function (Kolmogorov-smirnov)
    else if(type=="BHK"){
      LRI.statistic<-GoFBHK(x,a)
    }
    else stop(paste("unknown ", TYPE, "!"))
    return(list(statistic=LRI.statistic,lambda=e)) 
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="BH"){
    METHOD="Test of Baringhaus and Henze for the Exponential distribution"
  }else if(as.character(type)=="He"){
    METHOD="Test of Henze for the Exponential distribution"
  }else if(as.character(type)=="Kl"){
    METHOD="Test of Klar for the Exponential distribution"
  }
  else if(as.character(type)=="BHC"){
    METHOD="Test based on the integrated distribution function (Cramer-von mises) for the Exponential distribution"
  }
  else if(as.character(type)=="BHK"){
    METHOD="Test based on the integrated distribution function (Kolmogorov-smirnov) for the Exponential distribution"
  }
  stat <- LRI.statistic(x,type,a)
  statistic.obs <- stat$statistic
  estimate.obs <- stat$lambda
  fun<-function(y){
    fun <- LRI.statistic(y,type,a)
    return(fun$statistic)
  }
  sim.statistic <- GoFsim(nsim,n,fun) 
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  LRI.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                   method=paste(METHOD),
                   estimate = estimate.obs,  
                   data.name=DNAME)
  
  class(LRI.test) <-"htest"
  return(LRI.test)
}


CF.test<-function(x, type = "EP" , a = 1, nsim = 1000){
  #Family of the test statistics based on the caracteristic function 
  CF.statistic<-function(x, type, a){
    TYPE <- deparse(substitute(type))
    n = length(x)
    #EMV
    e = n/sum(x)    
    y=e*x
    l<- 1:n
    
    ### Test of Epps and Pulley
    if(type=="EP"){
      s <- mean(exp(-y)) -0.5 
      CF.statistic <- sqrt(48*n)*s  
    }
    ### Test of Henze and Meintanis 1
    else if(type=="W1"){
     s<-GoFHM(x,a)
     CF.statistic <-s[1]
    }
    ### Test of Henze and Meintanis 2
    else if(type=="W2"){

      s<-GoFHM(x,a)
      CF.statistic <-s[2]
    }
    ### Test of Meintanis and Iliopoulos 1
    else if(type=="T1"){
      s<-GoFMI(x,a)
      CF.statistic <-s[1]
    }
    
    ### Test of Meintanis and Iliopoulos 2
    else if(type=="T2"){
      s<-GoFMI(x,a)
      CF.statistic <-s[2]
    ### Test of Klar based on the integral of the distribution function
      } else if(type=="Kl"){
        ord_y=sort(y) 
        b=0
        c=0
        for (i in(1:(n-1))){ 
          b <- b+sum(ord_y[i]^2*ord_y[(i+1):n])  
          c <- c+sum(exp(-a*ord_y[i])*(a*(ord_y[(i+1):n]-ord_y[i])-2))         
        }
        s = 0
        if (a==0){
          s <- n/2-2*sum(exp(-y))-1/(3*n)*(t(n-l-1)%*%ord_y^3)+n^(-1)*b    
        } else { 
          s <- 2*(3*a+2)*n/((2+a)*(1+a)^2)-2*a^3/(1+a)^2*sum(exp(-(1+a)*y))-2/n*sum(exp(-a*y))+2/n*c 
        }     
        CF.statistic <- s
    }
    else stop(paste("unknown ", TYPE, "!"))
    return(list(statistic=CF.statistic,lambda=e)) 
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="EP"){
    METHOD="Test of Epps and Pulley for the Exponential distribution"
  } else if(as.character(type)=="W1"){
    if(a==0) stop(paste(TYPE, " is not defined for a=0 !"))
    METHOD="Test of Henze and Meintanis 1 for the Exponential distribution"
  } else if(as.character(type)=="W2"){
    if(a==0) stop(paste(TYPE, " is not defined for a=0 !"))
    METHOD="Test of Henze and Meintanis 2 for the Exponential distribution"
  } else if(as.character(type)=="T1"){
    if(a==0) stop(paste(TYPE, " is not defined for a=0 !"))
    METHOD="Test of Meintanis and Iliopoulos 1 for the Exponential distribution"
  } else if(as.character(type)=="T2"){
    if(a==0) stop(paste(TYPE, " is not defined for a=0 !"))
    METHOD="Test of Meintanis and Iliopoulos 2 for the Exponential distribution"
  }  else if(as.character(type)=="Kl"){
    METHOD="Test of Klar for the Exponential distribution"
  } 
  stat <- CF.statistic(x,type,a)
  statistic.obs <- stat$statistic
  estimate.obs <- stat$lambda
  fun<-function(y){
    fun <- CF.statistic(y,type,a)
    return(fun$statistic)
  }
  sim.statistic <- GoFsim(nsim,n,fun) 
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  CF.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                  method=paste(METHOD),
                  estimate = estimate.obs, 
                  data.name=DNAME)
  
  class(CF.test) <-"htest"
  return(CF.test)
}

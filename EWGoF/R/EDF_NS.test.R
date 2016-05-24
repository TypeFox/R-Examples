EDF_NS.test<-function(x,type="AD", nsim=2000){
  
  #Family of the test statistics based on the empirical distribution function 
  ECF.statistic<-function(x,type="AD"){
    TYPE <- deparse(substitute(type))
    n = length(x)
    #EMV
    e = n/sum(x)    
    #Computing the uniforms Ui
    u <- sort(x)
    l <- n+1-seq(1:n)
    #Normalized spacings
    E <- l*(u[1:n]-c(0,u[1:(n-1)]))   
    u <- 1-exp(-e*x)
    y <- sort(u)
    a <- 0
    l <- seq(1:n) 
    if(type=="CM"){
      ##Cramer
      a <- sum((y-(2*l-1)/(2*n))^2)
      a <- a + 1/(12*n)
      ECF.statistic <- a*(1 + 0.16/n)
    }
    ###Test statistic of Kolmogorov-Smirnov
    else if(type=="KS"){
      ##Kolmogorov Smirnov
      a <- max(max(l/n-y),max(y-(l-1)/n))
      a <- (a -0.2/n)*(sqrt(n)+0.26+0.5/sqrt(n))  
      ECF.statistic <- a
    }
    ###Test statistic of Anderson-Darling  
    else if(type=="AD"){
      ##Anderson-Darling
 
      a <- sum((2*l-1-2*n)*log(1-y)-(2*l-1)*log(y)) 
      a <- a/n -n
      ECF.statistic <- a*(1+0.6/n)
    ### Test statistic of Shapiro Wilk  
      } else if(type=="SW"){ 
        s2 <- sum(x^2)
        s <- sum(x)
        ECF.statistic <- s^2/(n*((n+1)*s2-s^2))
    ###Test statistic of Pathwardhan      
      } else if(type=="PA"){  
        ECF.statistic <- n*(n+1)*sum(E^2)/(sum(E))^2
    ### Test statistic of Gnedenko  
      }else if(type=="Gn"){
        r = floor(n/2)
        s = n-r
        ECF.statistic <- s*sum(E[1:r])/(r*sum(E[(r+1):n]))
      }else if(type=="G"){
        ECF.statistic <- sum(l[1:(n-1)]*E[2:n])/((n-1)*sum(E))
      }
    else stop(paste("unknown ", TYPE, "!"))
    return(list(statistic=ECF.statistic,lambda=e)) 
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="AD"){
    METHOD="Test of Anderson Darling for the Exponential distribution"
  } else if(as.character(type)=="CM"){
    METHOD="Test of Cramer von Mises for the Exponential distribution"
  } else if(as.character(type)=="KS"){
    METHOD="Test of Kolmogorov-Smirnov for the Exponential distribution"
  }  else if(as.character(type)=="SW"){
    METHOD="Test of Shapiro-Wlik for the Exponential distribution"
  } else if(as.character(type)=="PA"){
    METHOD="Test of Patwardhan for the Exponential distribution"
  } else if(as.character(type)=="Gn"){
    METHOD="Test of Gnedenko for the Exponential distribution"
  } else if(as.character(type)=="G"){
    METHOD="Test of Gini for the Exponential distribution"
  } 
  stat <- ECF.statistic(x,type)
  statistic.obs <- stat$statistic
  estimate.obs <- stat$lambda
  fun<-function(y){
    fun <- ECF.statistic(y,type)
    return(fun$statistic)
  }
  sim.statistic <- GoFsim(nsim,n,fun) 
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  pvalb <- 2*min(p_val,1-p_val)
  p_val <- switch(type,"SW"=pvalb,"G"=pvalb,p_val)
  EDF_NS.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                   method=paste(METHOD),
                   estimate = estimate.obs,  
                   data.name=DNAME)
  
  class(EDF_NS.test) <-"htest"
  return(EDF_NS.test)
}

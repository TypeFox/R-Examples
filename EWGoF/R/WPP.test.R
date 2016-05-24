

WPP.test <- function(x,type="SB",nsim=2000){
  ##Family of the test statistics based on the probability plot and shapiro-Wilk type tests
  WPP.statistic<-function(x,type="OK"){
    TYPE <- deparse(substitute(type))
    n <- length(x)
    lv <- log(x)
    y <- sort(lv)
    I <- 1:n
    #### Test statistic of Ozturk and Korukoglu
    if(type=="OK"){
      l <- I[1:(n-1)]
      Sig <- sum((2*c(l,n)-1-n)*y)/(.693147*(n-1))   
      w <- log((n+1)/(n-l+1))
      Wi <- c(w,n-sum(w))
      Wn <- w*(1+log(w))-1
      a <- 0.4228*n-sum(Wn)
      Wn <- c(Wn,a)
      b <- (0.6079*sum(Wn*y)-.2570*sum(Wi*y))
      stat <- b/Sig
      WPP.statistic <- (stat-1-.13/sqrt(n)+1.18/n)/(.49/sqrt(n)-.36/n)
      ####Test statistic of Shapiro Wilk
    }else if(type=="SB"){
      yb <- mean(y)
      S2 <- sum((y-yb)^2)
      l <- I[1:(n-1)]
      w <- log((n+1)/(n-l+1))
      Wi <- c(w,n-sum(w))
      Wn <- w*(1+log(w))-1
      a <- 0.4228*n-sum(Wn)
      Wn <- c(Wn,a)
      b <- (0.6079*sum(Wn*y)-.2570*sum(Wi*y))/n
      WPP.statistic <- n*b^2/S2
      ####Test statistic of Smith and Bain based on probability plot
    }else if(type=="RSB"){
      m <- I/(n+1)
      m <- log(-log(1-m))
      mb <- mean(m)
      xb <- mean(log(x))
      R <- (sum((log(x)-xb)*(m-mb)))^2
      R <- R/sum((log(x)-xb)^2)
      R <- R/sum((m-mb)^2)
      WPP.statistic <- n*(1-R)   
      ####Test statistic of Evans, Johnson and Green based on probability plot
    }else if(type=="REJG"){
      beta_shape <- MLEst(x)$beta
      m <- log(-(log(1-(I-.3175)/(n+0.365))))/beta_shape
      s <- (sum((y-mean(y))*m))^2/((sum((y-mean(y))^2))*sum((m-mean(m))^2))  
      WPP.statistic <- s^2
      ####Test statistic based on stabilized probability plot
    }else if(type=="SPP"){
      y=MLEst(x)$y
      r = 2/pi*asin(sqrt((I-0.5)/n))
      s = 2/pi*asin(sqrt(1-exp(-exp(y))))
      WPP.statistic <- max(abs(r-s))
      ### Smooth test statistic based on the skewness   
    }else if(type=="ST1"){
      x=sort(-y)
      s=sum((x-mean(x))^2)/n
      b1=sum(((x-mean(x))/sqrt(s))^3)/n
      V3=(b1-1.139547)/sqrt(20/n)
      WPP.statistic=V3^2
      ### Smooth test statistic based on the kurtosis
    }else if(type=="ST2"){
      x=sort(-y)
      s=sum((x-mean(x))^2)/n
      b1=sum(((x-mean(x))/sqrt(s))^3)/n
      b2=sum(((x-mean(x))/sqrt(s))^4)/n
      V4=(b2-7.55*b1+3.21)/sqrt(219.72/n)
      WPP.statistic=V4^2  
    }
    return(WPP.statistic)
  }
  DNAME <- deparse(substitute(x))
  TYPE <- deparse(substitute(type))
  n <- length(x)
  if(sum(x<0)){
    stop(paste("Data ", DNAME, " is not a positive sample"))}
  if(nsim<100){
    warning("small values of Monte-Carlo iterations may affect the value of the p-value")
  }
  if(as.character(type)=="OK"){
    METHOD="Test of Ozturk and Koruglu for the Weibull distribution"
  } else if(as.character(type)=="SB"){
    METHOD="Test of Shapiro and Bain for the Weibull distribution"
  }else if(as.character(type)=="RSB"){
    METHOD="Test of Smith and Bain for the Weibull distribution"
  }else if(as.character(type)=="REJG"){
    METHOD="Test of Evans, Johnson and Green for the Weibull distribution"
  }else if(as.character(type)=="SPP"){
    METHOD="Test based on the stabilized probability plot for the Weibull distribution"
  }else if(as.character(type)=="ST1"){
    METHOD="Smooth test statistic based on the skewness for the Weibull distribution"
  }else if(as.character(type)=="ST2"){
    METHOD="Smooth test statistic based on the kurtosis for the Weibull distribution"}
  else  stop(paste("The chosen method ",TYPE," is unknown"))
  
  MLE<-MLEst(x)
  estimate.obs <- c(MLE$eta,MLE$beta)
  statistic.obs <- WPP.statistic(x,type) 
  fun<-function(y){
    fun <- WPP.statistic(y,type)
    return(fun)
  }
  sim.statistic <- GoFsim(nsim,n,fun)
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  pvalb <- 2*min(p_val,1-p_val)
  p_val <- switch(type,"SB"=pvalb,"OK"=pvalb, p_val)
  WPP.test <- list(statistic =c(S=statistic.obs), p.value = p_val, 
                   method=paste(METHOD),
                   estimate = c(eta=estimate.obs[1],beta=estimate.obs[2]),  
                   data.name=DNAME)
  
  class(WPP.test) <-"htest"
  return(WPP.test)
}


WEDF.test<-function(x,type="AD",funEstimate="MLE",paramKL=2, nsim=2000){
 
  #Family of the test statistics based on the empirical distribution function 
  ECF.statistic<-function(x,type="AD",funEstimate="MLE",paramKL=2){
    TYPE <- deparse(substitute(type))
    funEst <- deparse(substitute(funEstimate))
    if(funEstimate=="MLE") {y.est<-MLEst(x)}
    else if (funEstimate=="ME"){y.est<-MEst(x)}
    else if (funEstimate=="LSE"){y.est<-LSEst(x)}
    #else if (funEstimate=="BLOM"){y.est<-BLOMEst(x)}
    else{stop(paste("unknown estimation method ", funEst, "!"))}
    z <- NULL  
    y<-y.est$y
    z <- exp(-exp(y)) 
    n<-length(x)
    l <- seq(length=n,from=1,to=n)
    z <- sort(z) 
    if(type=="CM"){
      ##Cramer
      res_CM <- sum((z-(2*l-1)/(2*n))^2)+ 1/(12*n)
      res_CM <- res_CM*(1+0.2/sqrt(n))
      ECF.statistic <- res_CM}
    else if(type=="W"){
      ##Cramer
      res_CM <- sum((z-(2*l-1)/(2*n))^2)+ 1/(12*n)
      ##Watson
      res_W <- res_CM-n*(mean(z)-0.5)^2
      res_W <- res_W*(1+0.2/sqrt(n))
      ECF.statistic <- res_W
    }
    ###Test statistic of Kolmogorov-Smirnoff
    else if(type=="KS"){
      ##Kolmogorov Smirnov
      res_KS <- 0
      res_KS <- max(max(l/n-z),max(z-(l-1)/n))
      ECF.statistic <- sqrt(n)*res_KS
    }
    ###Test statistic of Anderson-Darling  
    else if(type=="AD"){
      ##Anderson-Darling
      res_A <- 0
      d_z <- sort(z,decreasing=TRUE)
      res_A <- -n-1/n*sum((2*l-1)*(log(1-d_z)+log(z)))
      ECF.statistic <- res_A*(1+0.2/sqrt(n))
    }
    ###Test statistic of Liao-Shimokawa  
    else if(type=="LS"){ 
      ##Liao-Shimokawa
      res_LS <- 0
      res_A <- l
      res_A <- pmax((l/n-z),(z-(l-1)/n))
      res_A <- res_A/sqrt(z*(1-z))
      ECF.statistic <- sum(res_A)/sqrt(n)
    }
    ###Test statistic based on Kullback-Leibler information
    else if(type=="KL"){
      y <- sort(y)
      p<-0
      #Calcul de la statistique de test
      for (i in (1:n)){
        k = min(i+paramKL,n)
        t1 = y[k]
        l = max(i-paramKL,1)
        t2 = y[l]
        p = p -log(t1-t2)/n 
      }   
      ECF.statistic=p-mean(y)+mean(exp(y))-log(n/(2*paramKL))
      
      
    }
    else stop(paste("unknown ", TYPE, "!"))
    return(list(statistic=ECF.statistic,eta=y.est$eta,beta=y.est$beta)) 
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
  if(as.character(type)=="AD"){
    METHOD="Test of Anderson Darling for the Weibull distribution"
    family=1
  } else if(as.character(type)=="CM"){
    METHOD="Test of Cramer von Mises for the Weibull distribution"
    family=1
  } else if(as.character(type)=="KS"){
    METHOD="Test of Kolmogorov-Smirnoff for the Weibull distribution"
    family=1
  } else if(as.character(type)=="LS"){
    METHOD="Test of Liao and Shimokawa for the Weibull distribution"
    family=1
  } else if(as.character(type)=="W"){
    METHOD="Test of Watson for the Weibull distribution"
    family=1
  } else if(as.character(type)=="KL"){
    METHOD="Test based on Kullback-Leibler information for the Weibull distribution"
    #value of the parameter m of Kulback-Leibler information based test
    if(n>2 && n<=5){paramKL <- 2
    }else if(n>6 && n<=24){paramKL <- 3
    }else if(n>25 && n<=39){paramKL <- 4
    }else if(n>40 && n<=50){paramKL <- 5
    }else if(n>51 && n<=69){paramKL <- 6
    }else if(n>70 && n<=99){paramKL <- 7
    }else if(n>100 && n<=119){paramKL <- 8
    }else if(n>120 && n<=129){paramKL <- 9
    }else if(n>130 && n<=159){paramKL <- 10
    }else if(n>160 && n<=189){paramKL <- 12
    }else if(n>190 && n<=200){paramKL <- 13
    }else{paramKL<-14}
  }  
  EST <- deparse(substitute(funEstimate))
  if(as.character(funEstimate)=="MLE"){EST="using the MLEs "
  } else if(as.character(funEstimate)=="ME"){EST="using the MEs "
  } else if(as.character(funEstimate)=="LSE"){EST="using the LSEs "
  } else  stop(paste("The chosen estimation method ",EST," is unknown"))
  pramKL=switch(type,"KL"=paramKL,2)
  stat <- ECF.statistic(x,type,funEstimate,paramKL)
  statistic.obs <- stat$statistic
  estimate.obs <- c(stat$eta,stat$beta)
  fun<-function(y){
    fun <- ECF.statistic(y,type,funEstimate,paramKL)
    return(fun$statistic)
  }
  sim.statistic <- GoFsim(nsim,n,fun)
  p_val <- sum(sim.statistic>=statistic.obs)/nsim
  WEDF.test <- list(statistic =c(S=statistic.obs), p.value = p_val, method=paste(METHOD," ",EST),
                 estimate = c(eta=estimate.obs[1],beta=estimate.obs[2]),  
                 data.name=DNAME)
  
  class(WEDF.test) <-"htest"
  return(WEDF.test)
}

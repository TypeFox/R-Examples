testTVACD <- function(fitModel, K = 2, type = "total", robust = TRUE){  
  if(fitModel$model != "ACD" || fitModel$distribution != "exponential") stop("this test only works if the model estimated is a standard ACD model.")
  if(length(fitModel$durations$time) == 0) stop("this test requires 'fitModel' to have been estimated from a durationobject with the time of durations provided.")
  
  type <- match.arg(type, c("total", "intraday"))
  
  #initiates variables:
  p <- fitModel$order[1]
  q <- fitModel$order[2]  
  if(length(fitModel$durations$adjDur)){
    dur <- fitModel$durations$adjDur
  } else{
    dur <- fitModel$durations$durations    
  }
  maxpq <- max(p, q)
  N <- length(dur)
  k <- length(fitModel$mPara)
  mean <- mean(dur)
  dmudtheta<-matrix(nrow = length(dur), ncol = k) 
  dpsidtheta<-matrix(nrow = length(dur), ncol = (1+p+q)*K)
  mu <- fitModel$muHats
  beta <- fitModel$mPara[(2+p):k]
   
  if(type == "total"){
    time <- as.numeric(fitModel$durations$time - fitModel$durations$time[1])
    time <- time/max(time)    
  }else if(type == "intraday"){
    time <- fitModel$durations$time$hour * 3600 + fitModel$durations$time$min * 60 + fitModel$durations$time$sec
    time <- as.numeric(time - time[1])
    time <- time/max(time)
  }
    
  ######computes dmudtheta########
  zeros<-rep(0,maxpq)
  #for omega:
  dmudtheta[,1] <- c(zeros, stats::filter(rep(1,N-maxpq),beta,"r"))  
  #for alpha:
  for(j in 1:p){
    dmudtheta[,j+1] <- c(zeros, stats::filter(dur[(maxpq+1-j):(N-j)],beta,"r"))   
  }
  #for beta:  
  for(j in 1:q){
    dmudtheta[,j+p+1] <- c(zeros, stats::filter(mu[(maxpq+1-j):(N-j)],beta,"r"))
  }   
  
  ######computes dpsidtheta########
  #tl is needed for the partial derivatives:
  tl <- matrix(nrow = length(dur), ncol = K)
  for(l in 1:K){
    tl[ ,l] <- c(zeros, time[maxpq:(N-1)]^l)
  }
  #for the first K partial derivatives:
  for(l in 1:K){
    dpsidtheta[ ,l] <- c(zeros,stats::filter(tl[(maxpq+1):N, l],beta,"r"))
  }
    
  #for vec(dur1):
  for(l in 1:K){
    for(j in 1:p){
      dpsidtheta[, K*j+l] <- c(zeros, stats::filter(dur[(maxpq-j+1):(N-j)]*dpsidtheta[(maxpq+1):N, l],beta,"r"))
    }
  }  
  #for vec(dur2):
  for(l in 1:K){
    for(j in 1:q){
      dpsidtheta[, K*p+K*j+l] <- c(zeros, stats::filter(mu[(maxpq-j+1):(N-j)]*dpsidtheta[(maxpq+1):N, l],beta,"r"))
    }
  }
  
  ######computes a_i, b_i and c_i########
  a <- dmudtheta/mu
  b <- dpsidtheta/mu
  c <- dur/mu-1
  
  ######
  if(robust){ #"robust":
    
    regression1 <- stats::lm(b~a-1)
    cr <- c*regression1$residuals
    regression2 <- stats::lm(rep(1,N)~cr-1)
    SSR <- sum(regression2$residuals^2)      
    
    chi2 <- N - SSR
    pv <- 1 - stats::pchisq(chi2, (1+p+q)*K)    
    
  } else{ #the non robust version of the test:
    SSR0 <- sum(c^2)
    regression1 <- stats::lm(c ~ a + b - 1)
    SSR1 <- sum(regression1$residuals^2)      
    
    chi2 <- N * (SSR0 - SSR1) / SSR0
    pv <- 1 - stats::pchisq(chi2, (1+p+q)*K)      
  }
  
  df.out <- data.frame(c(chi2, (1+p+q)*K, pv))
  rownames(df.out) <- c("LM-stat: ", "Degrees of freedom: ", "P-value: ")
  colnames(df.out) <- " "
  
  if(robust) cat("\nM&T (2006) test of ACD against TVACD (robust version): \n")
  if(!robust) cat("\nM&T (2006) test of ACD against TVACD (non-robust version): \n")
  
  cat("\nType:",  type,"\n")
  print(format(df.out, digits = 3, scientific = F))
  
  testTVACD <- list(chi2 = chi2, pv = pv)
}
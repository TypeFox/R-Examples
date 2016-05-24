testSTACD<-function(fitModel, K = 2, robust = TRUE){
  
  if(fitModel$model != "ACD" || fitModel$distribution != "exponential") stop("this test only works if the model estimated is a standard ACD model.")
  
  #initiates variables:
  p <- fitModel$order[1]
  q <- fitModel$order[2]  
  if(length(fitModel$durations$adjDur)){
    dur <- fitModel$durations$adjDur
  } else{
    dur <- fitModel$durations$durations    
  }
  
  k <- sum(fitModel$order)+1
  maxpq <- max(p, q)
  N <- length(dur)
  mean <- mean(dur)
  dmudtheta<-matrix(nrow = N, ncol = k) 
  dpsidtheta<-matrix(nrow = N, ncol = 2*p*K)
  
  mu <- fitModel$muHats  
  beta <- fitModel$mPara[(2+p):k]
  
  ######computes dmudtheta########
  zeros<-rep(0,maxpq) 
  #for omega:
  dmudtheta[,1]<-c(zeros,stats::filter(rep(1,N-maxpq),beta,"r"))  
  #for alpha:
  for(j in 1:p){
    dmudtheta[,j+1]<-c(zeros,stats::filter(dur[(maxpq+1-j):(N-j)],beta,"r"))   
  }
  #for beta:  
  for(j in 1:q){
    dmudtheta[,j+p+1]<-c(zeros,stats::filter(mu[(maxpq+1-j):(N-j)],beta,"r"))
  }   
  
  ######computes dpsidtheta########
  lnX <- log(dur)
  for(l in 1:K){
    for(j in 1:p){
      dpsidtheta[,(j-1)*K+l]<-c(zeros,lnX[(maxpq+1-j):(N-j)]^l)
      dpsidtheta[,(j-1)*K+l]<-c(zeros,stats::filter(dpsidtheta[(maxpq+1):N,(j-1)*K+l],beta,"r"))
      dpsidtheta[,K*p+(j-1)*K+l]<-c(zeros,dur[(maxpq+1-j):(N-j)]*lnX[(maxpq+1-j):(N-j)]^l)
      dpsidtheta[,K*p+(j-1)*K+l]<-c(zeros,stats::filter(dpsidtheta[(maxpq+1):N,K*p+(j-1)*K+l],beta,"r"))
    }    
  }  
  
  ######computes a_i, b_i and c_i########
  a<-dmudtheta/mu
  b<-dpsidtheta/mu
  c<-dur/mu-1
  
  ######computes the LM-statistic########  
  if(robust){ #"robust":
    
    regression1 <- stats::lm(b~a-1)
    cr <- c*regression1$residuals
    regression2 <- stats::lm(rep(1,N)~cr-1)
    SSR <- sum(regression2$residuals^2)      
    
    chi2 <- N - SSR
    pv <- 1 - stats::pchisq(chi2, 2*p*K)    
    
  } else{ #the non robust version of the test:
    SSR0 <- sum(c^2)
    regression1 <- stats::lm(c ~ a + b - 1)
    SSR1 <- sum(regression1$residuals^2)      
    
    chi2 <- N * (SSR0 - SSR1) / SSR0
    pv <- 1 - stats::pchisq(chi2, 2*p*K)      
  }
  
  df.out <- data.frame(c(chi2, 2*p*K, pv))
  rownames(df.out) <- c("LM-stat: ", "Degrees of freedom: ", "P-value: ")
  colnames(df.out) <- " "
  
  if(robust) cat("\nM&T (2006) test of ACD against STACD (robust version): \n")
  if(!robust) cat("\nM&T (2006) test of ACD against STACD (nonrobust version): \n")
  print(format(df.out, digits = 3, scientific = F))
  
  testSTACD <- list(chi2 = chi2, pv = pv)
}
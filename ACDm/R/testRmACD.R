testRmACD <- function(fitModel, pStar = 2, robust = TRUE){
  if(fitModel$model != "ACD" || fitModel$distribution != "exponential") stop("this test only works if the model estimated is a standard ACD model.")
  
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
  dmudtheta <- matrix(nrow = N, ncol = k) 
  dpsidtheta <- matrix(nrow = N, ncol = pStar)
  mu <- fitModel$muHats
  beta <- fitModel$mPara[(2+p):k]
  
  ######comuptes dmudtheta########
  zeros<-rep(0, maxpq) 
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
  
  ######comuptes dpsidtheta########  
  zeros<-rep(0,pStar) 
  for(j in 1:pStar){
    dpsidtheta[ , j] <- c(zeros,dur[(pStar + 1 - j):(N-j)]/mu[(pStar + 1 - j):(N-j)])
  }  
  
  ######comuptes a_i, b_i och c_i########
  a <- dmudtheta/mu
  b <- dpsidtheta/mu
  c <- dur/mu-1
  
  ######comuptes the LM-statistic########
  
  if(robust){ #"robust":
    
    regression1 <- stats::lm(b~a-1)
    cr <- c*regression1$residuals
    regression2 <- stats::lm(rep(1,N)~cr-1)
    SSR <- sum(regression2$residuals^2)      
    
    chi2 <- N - SSR
    pv <- 1 - stats::pchisq(chi2, pStar)    
    
  } else{ #the non robust version of the test:
    SSR0 <- sum(c^2)
    regression1 <- stats::lm(c ~ a + b - 1)
    SSR1 <- sum(regression1$residuals^2)      
    
    chi2 <- N * (SSR0 - SSR1) / SSR0
    pv <- 1 - stats::pchisq(chi2, pStar)      
  }
  
  df.out <- data.frame(c(chi2, pStar, pv))
  rownames(df.out) <- c("LM-stat: ", "Degrees of freedom: ", "P-value: ")
  colnames(df.out) <- " "
  
  if(robust) cat("\nM&T (2006) test of no remaining ACD in residuals (robust version): \n")
  if(!robust) cat("\nM&T (2006) test of no remaining ACD in residuals (nonrobust version): \n")
  print(format(df.out, digits = 3, scientific = F))
  
  testRmACD <- list(chi2 = chi2, pv = pv)
}
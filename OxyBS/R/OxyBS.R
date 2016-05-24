# derivative of -log(beta pdf) wrt a parm
diffBeta1 <- function(x,a,b){
  digamma(a)-digamma(a+b)-log(x)
}

# derivative of -log(beta pdf) wrt b parm
diffBeta2 <- function(x,a,b){
  digamma(b)-digamma(a+b)-log(1-x)
}

# beta-likelihood function
likeOxBS <- function(theta, betaBS, betaOxBS, signalBS, signalOxBS){
  theta <- pmin(100,pmax(-100,theta))
  p <- exp(c(0,theta))
  p <- p/sum(p)
  a1  <- signalBS*(p[2]+p[3])
  b1  <- signalBS*(p[1])
  a2 <- signalOxBS*(p[2])
  b2 <- signalOxBS*(p[1]+p[3])  
  -(dbeta(betaBS,a1,b1,log=TRUE)+dbeta(betaOxBS,a2,b2,log=TRUE)) 
}

# Beta score function
scoreOxBS <- function(theta, betaBS, betaOxBS, signalBS, signalOxBS){
  theta <- pmin(100,pmax(-100,theta))
  p <- exp(c(0,theta))
  p <- p/sum(p)
  a1  <- signalBS*(p[2]+p[3])
  b1  <- signalBS*(p[1])
  a2 <- signalOxBS*(p[2])
  b2 <- signalOxBS*(p[1]+p[3])  

  dp <-   matrix(c(
     -p[1]*p[2], p[2]*(p[1]+p[3]), -p[2]*p[3],
     -p[1]*p[3], -p[2]*p[3], p[3]*(p[1]+p[2])
  ),3,2)

  ua1 <- ( diffBeta1(betaBS,a1,b1)*signalBS*(dp[2,1]+dp[3,1]) +
           diffBeta1(betaOxBS,a2,b2)*signalOxBS*(dp[2,1]) )
  ub1 <- ( diffBeta2(betaBS,a1,b1)*signalBS*(dp[1,1]) +
         diffBeta2(betaOxBS,a2,b2)*signalOxBS*(dp[1,1]+dp[3,1]) )

  ua2 <- ( diffBeta1(betaBS,a1,b1)*signalBS*(dp[2,2]+dp[3,2]) +
           diffBeta1(betaOxBS,a2,b2)*signalOxBS*(dp[2,2]) )
  ub2 <- ( diffBeta2(betaBS,a1,b1)*signalBS*(dp[1,2])+
           diffBeta2(betaOxBS,a2,b2)*signalOxBS*(dp[1,2]+dp[3,2]) )

  c(ua1+ub1,ua2+ub2)

}

# Fit one value
fitOneOxBS <- function(betaBS, betaOxBS, signalBS, signalOxBS, eps=1E-5){
  est5mC <- max(min(betaOxBS,1-eps),eps)
  estTotMeth <- max(min(betaBS,1-eps),eps)
  est5hmC <- max(estTotMeth-est5mC,eps)
  
  #theta <- log(c(est5mC,est5hmC))-log(1-estTotMeth)  # Prior to Nov. 30, 2015
  theta <- log(c(est5mC,est5hmC))-log(1-est5mC-est5hmC)  # Fix after Nov. 30, 2015
  opt <- try( optim(theta, likeOxBS, method="BFGS", gr=scoreOxBS, 
    betaBS=betaBS, betaOxBS=betaOxBS, 
    signalBS=signalBS, signalOxBS=signalOxBS) )

  if(inherits(opt,"try-error")){
    #print(c(est5mC,est5hmC,estTotMeth))
    out <- c(1-estTotMeth,est5mC,est5hmC)
  }
  else out <- exp(c(0,opt$par))

  #names(out) <- c("C","5mC","5hmC")
  out/sum(out)
}

# Fit multiple values
fitOxBS <- function(betaBS, betaOxBS, signalBS, signalOxBS, eps=1E-5){
  n <- length(betaBS)
  out <- matrix(NA,n,3)
  for(i in 1:n){
    out[i,] <- fitOneOxBS(betaBS[i],betaOxBS[i], signalBS[i], signalOxBS[i],eps=eps)
    if(i %% 5000==0) cat(i,"\n")
  }
  colnames(out) <- c("C","5mC","5hmC")
  rownames(out) <- names(betaBS)
  out
}


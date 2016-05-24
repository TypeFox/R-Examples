###################################################################################
### Purpose: Find the MLE of log marginalized likelihood function from
###          Sarmanov beta distribution or independent beta distribution
### Input:   data(y1,n1,y2,n2), model
### Output:  a list, containning pvalue,chi2,MLE and 
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
MLE.function <- function(y1=y1,n1=n1,y2=y2,n2=n2,model=model) {
  init.val <- rep(0, 5)
  fit1 <- initial.val.gen(y1, n1)                                                 
  fit2 <- initial.val.gen(y2, n2)
  init.val[1] <- log(fit1$a);                                                     
  init.val[2] <- log(fit1$b)   
  init.val[3] <- log(fit2$a);
  init.val[4] <- log(fit2$b)

  MLE.sar.log <- optim(init.val, myLik.sar.log, method = "L-BFGS-B",
                       lower=rep(-20,5), upper=rep(20,5),
                       control = list(fnscale=-1,maxit=1000),
                       hessian=T, mydat=list(y1=y1,n1=n1,y2=y2,n2=n2))
  MLE.inde.log <- optim(init.val[1:4], myLik.indep.log, method = "L-BFGS-B",
                        lower=rep(-20,4), upper=rep(20,4),
                        control = list(fnscale=-1,maxit=1000),
                        hessian = T, mydat=list(y1=y1,n1=n1,y2=y2,n2=n2))

  if (model=="Sarmanov") {
    mypar <- par.cal(MLE.sar.log$par); ##tranfer back to original scale
    rho <- mypar[5];
    hessian.log <- MLE.sar.log$hessian
    colnames(hessian.log) <- c("loga1","logb1","loga2","logb2","eta")
    rownames(hessian.log) <- c("loga1","logb1","loga2","logb2","eta")
  }
  if (model=="Independent") {
    mypar <- par.cal(MLE.inde.log$par)
    rho <- 0
    hessian.log <- MLE.inde.log$hessian
    colnames(hessian.log) <- c("loga1","logb1","loga2","logb2")
    rownames(hessian.log) <- c("loga1","logb1","loga2","logb2")
  }                                          
  a1 <- mypar[1]; b1 <- mypar[2]
  a2 <- mypar[3]; b2 <- mypar[4]
  prior.MLE <- c(a1, b1, a2, b2, rho)
  chi2 <- (-2*(MLE.inde.log$value - MLE.sar.log$value))
  p.value <- round(pchisq(-2*(MLE.inde.log$value - MLE.sar.log$value),
                          df=1, lower.tail=F),2)
 return(list(chi2=chi2, pvalue=p.value,MLE=prior.MLE,hessian.log=hessian.log))                                                                           
}

logit <- function(p) log(p/(1-p))
expit <- function(x) exp(x)/(1+exp(x))

##################################################################################################
### Purpose: Compute the log marginalized likelihood function from Sarmanov beta distribution 
### Input: 1)mypar: log(a1,b1,a2,b2), eta
###        2)mydata: n1,y1,y2,n2
### Output:  loglikilhood
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################################
myLik.sar.log <- function(mypar, mydat) { 
  a1.temp <- mypar[1]; b1.temp <- mypar[2]
  a2.temp <- mypar[3]; b2.temp <- mypar[4]
  eta <- mypar[5]

  a1 <- exp(a1.temp); b1 <- exp(b1.temp)
  a2 <- exp(a2.temp); b2 <- exp(b2.temp)

  # mu1, mu2: means
  mu1 <- a1/(a1+b1); mu2 <- a2/(a2+b2)  
  # delta1, delta2: standard deviations
  delta1 <- sqrt(mu1*(1-mu1)/(a1+b1+1))
  delta2 <- sqrt(mu2*(1-mu2)/(a2+b2+1))
  
  cc <- sqrt(a1*a2*b1*b2)/sqrt((a1+b1+1)*(a2+b2+1))
  upper.bound <- cc/max(a1*b2, a2*b1)
  lower.bound <- -cc/max(a1*a2, b1*b2)
  
  rho <- (upper.bound-lower.bound)*expit(eta) + lower.bound

  temp1 <- (lgamma(a1+mydat$y1) + lgamma(b1+mydat$n1-mydat$y1)
            + lgamma(a2+mydat$y2) + lgamma(b2+mydat$n2-mydat$y2)
            + lgamma(a1+b1) + lgamma(a2+b2))
  temp2 <- (lgamma(a1) + lgamma(b1) + lgamma(a2) + lgamma(b2)
            + lgamma(a1+b1+mydat$n1) + lgamma(a2+b2+mydat$n2))
  temp3 <- (log(1+rho/delta1/delta2
                        *(mydat$y1-mydat$n1*mu1)
                        *(mydat$y2-mydat$n2*mu2)
                        /(a1+b1+mydat$n1)/(a2+b2+mydat$n2)))
  myLogLik <- sum(temp1) - sum(temp2) + sum(temp3)

  return(myLogLik)  
}

################################################################################################# 
### Purpose: Compute the log marginalized likelihood function from independent beta distribution
### Input: 1)mypar: log(a1,b1,a2,b2), eta
###        2)mydata: n1,y1,y2,n2
### Output:  loglikilhood
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
##################################################################################################
myLik.indep.log <- function(mypar, mydat) {
  a1.temp <- mypar[1]; b1.temp <- mypar[2]
  a2.temp <- mypar[3]; b2.temp <- mypar[4]

  a1 <- exp(a1.temp); b1 <- exp(b1.temp)
  a2 <- exp(a2.temp); b2 <- exp(b2.temp)

  temp1 <- (lgamma(a1+mydat$y1) + lgamma(b1+mydat$n1-mydat$y1)
            + lgamma(a2+mydat$y2) + lgamma(b2+mydat$n2-mydat$y2)
            + lgamma(a1+b1) + lgamma(a2+b2))
  temp2 <- (lgamma(a1) + lgamma(b1) + lgamma(a2) + lgamma(b2)
            + lgamma(a1+b1+mydat$n1) + lgamma(a2+b2+mydat$n2))
  
  myLogLik <- sum(temp1 - temp2)

  return(myLogLik)  
}

########################################################################################################
### Purpose: This function fit beta-binomial model to generate initial values for hyperparameters a & b
### Input: y and n. Both are data vector. y is the number of events. n is the number of experiments
### Output: initial value of hyperparemeters a1,b1,a2,b2 for opzimazation
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
########################################################################################################
initial.val.gen <- function(y, n) {
  BBfit <- aod::betabin(cbind(y, n-y)~1, ~1, data=data.frame(y=y,n=n))
  a.ini <- expit(as.numeric(BBfit@param[1]))*(1/as.numeric(BBfit@param[2])-1)
  b.ini <- (1/as.numeric(BBfit@param[2])-1)*(1-expit(as.numeric(BBfit@param[1])))
  return(list(a=a.ini, b=b.ini))
}

###################################################################################################
### Purpose: This function calculates the inverse of a symmetric matrix (e.g. Hessian/Information)
###           It is standardized first to avoid computational problem
### Input: symmetric matrix: Hessian/Information : parts of the results of optim
### Output: inverse matrix or warming
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
inverse.matrix.func <- function(Mat){
    ## Mat: a square symmetric matrix with positive diagonal elements
    ## checking Mat to be correctly defined
    ## (i.e. no NA/Infinite, symmetric, no negative or zero diagonal element)
    ## if not, the inverse is a matrix of the same dimension of Mat with NAs

  Mat.inverse <- matrix(NA, nrow=nrow(Mat), ncol=ncol(Mat))

  not.defined.logic <- (((sum(is.na(Mat))+ sum(is.infinite(Mat)))>0)|(mean(t(Mat)-Mat)>1e-4)|(sum(diag(Mat)<=0)>0))
  if(!not.defined.logic){
    Mat.inverse <- matrix(NA, nrow=nrow(Mat), ncol=ncol(Mat))
    Mat.stand <- diag(1/sqrt(diag(Mat)))%*%Mat%*%diag(1/sqrt(diag(Mat)))
    Mat.rank <- qr(Mat.stand)$rank
    if(Mat.rank==nrow(Mat)){
      Mat.inverse <- diag(1/sqrt(diag(Mat)))%*%solve(Mat.stand)%*%diag(1/sqrt(diag(Mat)))
    }
  }
  Mat.inverse
}

##################################################################################################
### Purpose:  Transfer the paremeter estimates in the transformed scales into the original scale
### input: parameters in the transformed scales: log(a1,b1,a2,b2), eta
### Output: parameters in the original scale
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################################
par.cal <- function(mypar) {
  a1 <- exp(mypar[1]); b1 <- exp(mypar[2])
  a2 <- exp(mypar[3]); b2 <- exp(mypar[4])
  eta <- mypar[5]
  cc <- sqrt(a1*a2*b1*b2)/sqrt((a1+b1+1)*(a2+b2+1))
  upper.bound <- cc/max(a1*b2, a2*b1)
  lower.bound <- -cc/max(a1*a2, b1*b2)
  rho <- (upper.bound-lower.bound)*expit(eta) + lower.bound
  return(c(a1,b1,a2,b2,rho))
}

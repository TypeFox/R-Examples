###################################################################################
### Purpose: Transform parameter from original scale to logit scale
### Input:   original scale
### Output:  logit scale
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
logit <- function(p) log(p/(1-p))

###################################################################################
### Purpose: Transform parameter from logit scale to original scale
### Input:   logit scale
### Output:  original scale
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
expit <- function(x) exp(x)/(1+exp(x))




###################################################################################
### Purpose: Compute confidence interval from samples
### Input:   a vector contains samples, significance level
### Output:  list, contraining confidence interval
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################
minlength.CI <- function(x, alpha, left=-1E10, right=1E10){                     
  n <- length(x)                                                                 
  sortx <- sort(c(x, left, right))
  disx <- sortx[2:(n+2)]-sortx[1:(n+1)]
  xn <- ceiling((n+1)*alpha)
  xnsum <- cbind(1:xn, cumsum(disx[1:xn])+sum(disx[(n-xn+2):(n+1)])-cumsum(c(0, disx[(n-xn+2):n])))
  nleft <- mean(xnsum[xnsum[, 2]==max(xnsum[, 2]), 1])
  CItemp <- quantile(sortx, probs=c(alpha*(nleft-1)/(xn-1), 1-alpha+alpha*(nleft-1)/(xn-1)))
  return(CItemp)
}
 
###################################################################################
### Purpose: compute kth moment of OR for independent model
### Input:   k, hyperparemeters(a1, b1, a2, b2)
### Output:  the kth moment
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
################################################################################### 
moment.OR.inde <- function(k, alpha1, beta1, alpha2, beta2) {
  mylog <- ((lgamma(alpha1-k)+lgamma(beta1+k)+lgamma(alpha2+k)+lgamma(beta2-k))
           -(lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)+lgamma(beta2)))
  return(exp(mylog))
}
                                                                  
###################################################################################
### Purpose: compute kth moment of OR for Sarmanov model
### Input:   k, hyperparemeters(a1,b1,a2,b2)
### Output:  the kth moment
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
################################################################################### 
moment.OR.sar <- function(k, y1, n1, y2, n2, a1, b1, a2, b2, rho) {
  myOmega <- omega.computation(y1, n1, y2, n2, a1, b1, a2, b2, rho)
  alpha1 <- y1+a1; beta1 <- n1-y1+b1
  alpha2 <- y2+a2; beta2 <- n2-y2+b2  
  return(myOmega$omega1*moment.OR.inde(k, alpha1, beta1, alpha2, beta2)
        +myOmega$omega2*moment.OR.inde(k, alpha1+1, beta1, alpha2, beta2)
        +myOmega$omega3*moment.OR.inde(k, alpha1, beta1, alpha2+1, beta2)
        +myOmega$omega4*moment.OR.inde(k, alpha1+1, beta1, alpha2+1, beta2))   
}

###################################################################################
### Purpose: compute kth moment of RR for independent model
### Input:   k, hyperparemeters(a1,b1,a2,b2)
### Output:  the kth moment
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
################################################################################### 
moment.RR.inde <- function(k, alpha1, beta1, alpha2, beta2) {
  mylog <- ((lgamma(alpha1-k)+lgamma(beta1)+lgamma(alpha2+k)+lgamma(beta2))
           -(lgamma(alpha1)+lgamma(beta1)+lgamma(alpha2)+lgamma(beta2)))
  return(exp(mylog))
                                                            }                                  
###################################################################################
### Purpose: compute kth moment of RR for Sarmanov model
### Input:   k, hyperparemeters(a1,b1,a2,b2)
### Output:  the kth moment
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
###################################################################################                                                             
moment.RR.sar <- function(k, y1, n1, y2, n2, a1, b1, a2, b2, rho) {
  myOmega <- omega.computation(y1, n1, y2, n2, a1, b1, a2, b2, rho)
  alpha1 <- y1+a1; beta1 <- n1-y1+b1
  alpha2 <- y2+a2; beta2 <- n2-y2+b2  
  return(myOmega$omega1*moment.RR.inde(k, alpha1, beta1, alpha2, beta2)
        +myOmega$omega2*moment.RR.inde(k, alpha1+1, beta1, alpha2, beta2)
        +myOmega$omega3*moment.RR.inde(k, alpha1, beta1, alpha2+1, beta2)
        +myOmega$omega4*moment.RR.inde(k, alpha1+1, beta1, alpha2+1, beta2))   
}                    

################################################################################################
### Purpose: compute kth moment of OR/RR for independent/Sarmanov model: Wrapper function
### Input:   k, hyperparemeters(a1,b1,a2,b2,rho),measure,model
### Output:  the kth moment
### Note:    Implemented by "moment.OR.inde","moment.OR.sar","moment.RR.inde","moment.RR.sar "
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
################################################################################################ 
moment.function <- function(measure=measure, model=model, k=k, y1=y1, n1=n2, y2=y2,
                            n2=n2, a1=a1, b1=b1, a2=a2, b2=b2, rho=rho){
  alpha1 <- y1+a1; beta1 <- n1-y1+b1
  alpha2 <- y2+a2; beta2 <- n2-y2+b2
  if(model=="Independent") {
    if (measure=="OR")  result <- moment.OR.inde(k, alpha1, beta1, alpha2, beta2)
    if (measure=="RR")  result <- moment.RR.inde(k, alpha1, beta1, alpha2, beta2)
  }
  if (model=="Sarmanov") {
    if (measure=="OR")  result <- moment.OR.sar(k, y1, n1, y2, n2, a1, b1, a2, b2, rho)
    if (measure=="RR")  result <- moment.RR.sar(k, y1, n1, y2, n2, a1, b1, a2, b2, rho)
  }
  return(result)
}
 

 
 

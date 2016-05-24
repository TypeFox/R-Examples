##############################################################################################
### Purpose: Wrapper function for computation of overall comparative measure
### Input:   MLE(a1,a2,b1,b2,rho), hessian matrix(in transformed scale), comparative measure,
###          model, alpha
### Output:  list, containning point estimate and CI, in original scale
### Note:    This function uses functions "OR.sar.log", "RR.sar.log", "RD.sar"
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
##############################################################################################
overall <- function(MLE, hessian.log, measure, model, alpha=alpha) {
  a1 <- MLE[1]; b1 <- MLE[2]
  a2 <- MLE[3]; b2 <- MLE[4]
  rho <- MLE[5]
  if (measure=="OR") ## hessian is in log scale
    out <- OR.sar.log(a1, b1, a2, b2, hessian.log,alpha=alpha)
  if (measure=="RR")
    out <- RR.sar.log(a1, b1, a2, b2, hessian.log,alpha=alpha)
  if (measure=="RD")
    out <- RD.sar(a1, b1, a2, b2, hessian.log,alpha=alpha) #point estimate of overal RD # wald CI of
  return(out)                                                   
}

#############################################################################################
### Purpose: Compute point estiamte and wald CI of the overall OR
### Input:   MLE(a1,a2,b1,b2) in log scale, hessian matrix(in log scale), significance level
### Output:  list contains point estimate and wald CI in original scale
### Note:    This function is used in wrapper function "overall". Variance is computated
###          via Delta method.
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
############################################################################################# 
OR.sar.log <- function(a1, b1, a2, b2, hessian.log, alpha=alpha) {
  myOR.log <- log((a2/b2)/(a1/b1))
  myVar.log <- inverse.matrix.func(-hessian.log)  
  
 if (nrow(hessian.log)==4) myD <- matrix(c(-1, 1, 1, -1), nrow=1)  else
  myD <- matrix(c(-1, 1, 1, -1, 0), nrow=1)
#
  myVar.log <- as.numeric(myD %*% myVar.log %*% t(myD))
  myOR.log.left <- myOR.log-qnorm(1-alpha/2)*sqrt(myVar.log)
  myOR.log.right <- myOR.log +qnorm(1-alpha/2)*sqrt(myVar.log)
  return(list(overall=exp(myOR.log),
              CI=exp(c(myOR.log.left,myOR.log.right))))
}                                                     

#############################################################################################
### Purpose: Compute point estiamte and wald CI of the overall RR
### Input:   MLE(a1,a2,b1,b2) in log scale, hessian matrix(in log scale), significance level
### Output:  list contains point estimate and wald CI in original scale
### Note:    This function is used in wrapper function "overall". Variance is computated
###          via Delta method.
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
#############################################################################################
RR.sar.log <- function(a1, b1, a2, b2, hessian.log,alpha=alpha) {
  myRR.log <- log((a2/(a2+b2))/(a1/(a1+b1)))
  myVar <- inverse.matrix.func(-hessian.log)
  if (nrow(hessian.log)==4)
    myD <- matrix(c(-b1/(a1+b1), b1/(a1+b1), b2/(a2+b2), -b2/(a2+b2)), nrow=1)
  if (nrow(hessian.log)==5)    
  myD <- matrix(c(-b1/(a1+b1), b1/(a1+b1), b2/(a2+b2), -b2/(a2+b2), 0), nrow=1)  

  myRR.log.Var <- as.numeric(myD %*% myVar %*% t(myD))
  myRR.log.left<- myRR.log-qnorm(1-alpha/2)*sqrt(myRR.log.Var)                       
  myRR.log.right <- myRR.log+qnorm(1-alpha/2)*sqrt(myRR.log.Var)                      
  return(list(overall=exp(myRR.log), CI=exp(c(myRR.log.left,myRR.log.right))))
}

#############################################################################################
### Purpose: Compute point estiamte and wald CI of the overall RD
### Input:   MLE(a1,a2,b1,b2) in log scale, hessian matrix(in log scale), significance level
### Output:  list contains point estimate and wald CI in original scale
### Note:    This function is used in wrapper function "overall". Variance is computated
###          via Delta method.
### Author:  Sheng Luo, Yong Chen, Xiao Su, Haitao Chu
### Data:    7/13/2012
#############################################################################################
RD.sar <- function(a1, b1, a2, b2, hessian.log,alpha=alpha) {
  myRD <- a2/(a2+b2)-a1/(a1+b1) 
  myVar <- inverse.matrix.func(-hessian.log)
  if (nrow(hessian.log)==4)
    myD <- matrix(c(-a1*b1/(a1+b1)^2, a1*b1/(a1+b1)^2, a2*b2/(a2+b2)^2,
                    -a2*b2/(a2+b2)^2), nrow=1)
  else myD <- matrix(c(-a1*b1/(a1+b1)^2, a1*b1/(a1+b1)^2, a2*b2/(a2+b2)^2,
                       -a2*b2/(a2+b2)^2, 0), nrow=1)
  myRD.Var <- as.numeric(myD %*% myVar %*% t(myD))
  myRD.left <- max(myRD-qnorm(1-alpha/2)*sqrt(myRD.Var), -1)
  myRD.right<- min(myRD+qnorm(1-alpha/2)*sqrt(myRD.Var), 1)
  # output of RD is in log scale
  return(list(overall=myRD, CI=c(myRD.left,myRD.right)))  
}

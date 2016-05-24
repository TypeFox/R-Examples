
##########################################################################
###########  3 STATES
##########################################################################


##########################################################################
###########  Tak K Mak
##########################################################################
##### Second derivative of the cost function (3 states): just a check
##### arguments:  x=thres3 object
##### value:      the second derivative
##########################################################################
secondDer3 <- function(x){
  # error handling
  if (class(x) != "thres3"){
    stop("'x' must be of class 'thres3'")
  }  
  if (x$T$dist1=="norm" & x$T$dist2=="norm" & x$T$dist3=="norm"){
    k1 <- x$T$k1
    k2 <- x$T$k2
    k3 <- x$T$k3
    rho <- x$T$prev
    costs <- x$T$costs
    Thr <- c(x$T$thres1, x$T$thres2)
    # sort distributions 
    mean <- c(mean(k1), mean(k2), mean(k3))
    ord_m <- order(mean)
    if (!all(ord_m==1:3)){
      list.ks.ord <- list(k1=k1, k2=k2, k3=k3)[ord_m]
      k1 <- list.ks.ord[[1]]
      k2 <- list.ks.ord[[2]]
      k3 <- list.ks.ord[[3]]
      rho <- rho[ord_m]
      costs <- costs[ord_m, ord_m]
    }
    # end of sorting    
    par1.1 <- mean(k1)
    par2.1 <- mean(k2)
    par3.1 <- mean(k3)
    par1.2 <- sd(k1)
    par2.2 <- sd(k2)
    par3.2 <- sd(k3)
    n1 <- length(k1)
    n2 <- length(k2)
    n3 <- length(k3)
    rho1 <- rho[1]
    rho2 <- rho[2]
    rho3 <- rho[3]
    c11 <- costs[1,1]
    c12 <- costs[1,2]
    c13 <- costs[1,3]
    c21 <- costs[2,1]
    c22 <- costs[2,2]
    c23 <- costs[2,3]
    c31 <- costs[3,1]
    c32 <- costs[3,2]
    c33 <- costs[3,3]  
    der1 <- (n1+n2+n3)/sqrt(2*pi)*( rho1*(c12-c11)*((Thr[1]-par1.1)/par1.2^3)*exp(-(Thr[1]-par1.1)^2/(2*par1.2^2))+
                                    rho2*(c22-c21)*((Thr[1]-par2.1)/par2.2^3)*exp(-(Thr[1]-par2.1)^2/(2*par2.2^2))+
                                    rho3*(c32-c31)*((Thr[1]-par3.1)/par3.2^3)*exp(-(Thr[1]-par3.1)^2/(2*par3.2^2)) )
    der2 <- (n1+n2+n3)/sqrt(2*pi)*( rho1*(c13-c12)*((Thr[2]-par1.1)/par1.2^3)*exp(-(Thr[2]-par1.1)^2/(2*par1.2^2))+
                                    rho2*(c23-c22)*((Thr[2]-par2.1)/par2.2^3)*exp(-(Thr[2]-par2.1)^2/(2*par2.2^2))+
                                    rho3*(c33-c32)*((Thr[2]-par3.1)/par3.2^3)*exp(-(Thr[2]-par3.1)^2/(2*par3.2^2)) )
  }else{
    dist1 <- x$T$dist1
    dist2 <- x$T$dist2
    dist3 <- x$T$dist3
    par1.1 <- x$T$pars1[1]
    par1.2 <- x$T$pars1[2]
    par2.1 <- x$T$pars2[1]
    par2.2 <- x$T$pars2[2]
    par3.1 <- x$T$pars3[1]
    par3.2 <- x$T$pars3[2]
    rho <- x$T$prev
    costs <- x$T$costs
    # sort distributions
    median1 <- quant(dist1)(0.5, par1.1, par1.2)
    median2 <- quant(dist2)(0.5, par2.1, par2.2)
    median3 <- quant(dist3)(0.5, par3.1, par3.2)  
    median <- c(median1, median2, median3)
    ord_m <- order(median)
    if (!all(ord_m==1:3)){
      pars <- matrix(c(par1.1, par1.2, par2.1, par2.2, par3.1, par3.2), byrow=T, ncol=2)
      pars <- pars[ord_m,]
      rho <- rho[ord_m]
      costs <- costs[ord_m,ord_m]
      par1.1 <- pars[1, 1]
      par1.2 <- pars[1, 2]
      par2.1 <- pars[2, 1]
      par2.2 <- pars[2, 2]
      par3.1 <- pars[3, 1]
      par3.2 <- pars[3, 2]
      dist <- c(dist1, dist2, dist3)[ord_m]
      dist1 <- dist[1]
      dist2 <- dist[2]
      dist3 <- dist[3]
    }
    # end of sorting
    rho1 <- rho[1]
    rho2 <- rho[2]
    rho3 <- rho[3]
    c11 <- costs[1, 1]
    c12 <- costs[1, 2]
    c13 <- costs[1, 3]
    c21 <- costs[2, 1]
    c22 <- costs[2, 2]
    c23 <- costs[2, 3]
    c31 <- costs[3, 1]
    c32 <- costs[3, 2]
    c33 <- costs[3, 3] 
    eq1 <- function(y){
      ret <- (rho1*(c11-c12)*dens(dist1)(y,par1.1,par1.2)+rho2*(c21-c22)*dens(dist2)(y,par2.1,par2.2)+rho3*(c31-c32)*dens(dist3)(y,par3.1,par3.2))
      return(ret)
    }
    eq2 <- function(y){
      ret <- (rho1*(c12-c13)*dens(dist1)(y,par1.1,par1.2)+rho2*(c22-c23)*dens(dist2)(y,par2.1,par2.2)+rho3*(c32-c33)*dens(dist3)(y,par3.1,par3.2))
      return(ret)
    }
    der1 <- grad(eq1, x$T$thres1)
    der2 <- grad(eq2, x$T$thres2)
  }
  der <- c(der1, der2)
  names(der) <- c("Value for thres1", "Value for thres2")
  return(der)
}

##########################################################################
##### Second derivative of the cost function (3 states): just a check
##### arguments:  k1=vector containing the class 1 sample values 
#####  	   	      k2= vector containing the class 2 sample values
#####		   	      k3= vector containing the class 3 sample values 
#####	 	   	      rho= vector of prevalences
#####		  	      costs=cost matrix
#####		  	      Thr=at which threshold values will the second derivative be evaluated
#####             na.rm=a logical value indicating whether NA values in k1, k2 and k3 should
#####             be stripped before the computation proceeds. Default, F.
##### value:      the second derivative function
##########################################################################
secondDer3aux <- function(k1, k2, k3, rho, Thr, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE), na.rm=FALSE){
  # error handling
  if (sum(rho > 0 & rho < 1) != 3){
    stop("The prevalences must be in (0,1)")
  }
  if (sum(rho) != 1){
    stop("The sum of the prevalences must be 1")
  }
  if (!is.matrix(costs)){
    stop("'costs' must be a matrix")
  }
  if (dim(costs)[1] != 3 | dim(costs)[2] != 3){
    stop("'costs' must be a 3x3 matrix")
  }
  if (!is.numeric(k1) | !is.numeric(k2) | !is.numeric(k3)){
    stop("'k1', 'k2' and 'k3' must be numeric vectors")
  }
  # NAs handling
  if (na.rm){
    k1 <- k1[!is.na(k1)]
    k2 <- k2[!is.na(k2)]
    k3 <- k2[!is.na(k3)]
  }
  par1.1 <- mean(k1)
  par2.1 <- mean(k2)
  par3.1 <- mean(k3)
  par1.2 <- sd(k1)
  par2.2 <- sd(k2)
  par3.2 <- sd(k3)
  n1 <- length(k1)
  n2 <- length(k2)
  n3 <- length(k3)
  rho1 <- rho[1]
  rho2 <- rho[2]
  rho3 <- rho[3]
  c11 <- costs[1, 1]
  c12 <- costs[1, 2]
  c13 <- costs[1, 3]
  c21 <- costs[2, 1]
  c22 <- costs[2, 2]
  c23 <- costs[2, 3]
  c31 <- costs[3, 1]
  c32 <- costs[3, 2]
  c33 <- costs[3, 3]
  
  der1 <- (n1+n2+n3)/sqrt(2*pi)*( rho1*(c12-c11)*((Thr[1]-par1.1)/par1.2^3)*exp(-(Thr[1]-par1.1)^2/(2*par1.2^2))+
                                    rho2*(c22-c21)*((Thr[1]-par2.1)/par2.2^3)*exp(-(Thr[1]-par2.1)^2/(2*par2.2^2))+
                                    rho3*(c32-c31)*((Thr[1]-par3.1)/par3.2^3)*exp(-(Thr[1]-par3.1)^2/(2*par3.2^2)) )
  der2 <- (n1+n2+n3)/sqrt(2*pi)*( rho1*(c13-c12)*((Thr[2]-par1.1)/par1.2^3)*exp(-(Thr[2]-par1.1)^2/(2*par1.2^2))+
                                    rho2*(c23-c22)*((Thr[2]-par2.1)/par2.2^3)*exp(-(Thr[2]-par2.1)^2/(2*par2.2^2))+
                                    rho3*(c33-c32)*((Thr[2]-par3.1)/par3.2^3)*exp(-(Thr[2]-par3.1)^2/(2*par3.2^2)) )
  der <- c(der1, der2)
  names(der) <- c("Value for Thr[1]","Value for Thr[2]")
  return(der)
}


##########################################################################
###########  V=Delta Variances of the derivatives of the cost function sigma^2/2*n
#####		  	costs=cost matrix
##########################################################################
VDel3 <- function(k1, k2, k3, rho, Thr, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE)){
  n1 <- length(k1)
  n2 <- length(k2)
  n3 <- length(k3)
  c11 <- costs[1, 1]
  c12 <- costs[1, 2]
  c13 <- costs[1, 3]
  c21 <- costs[2, 1]
  c22 <- costs[2, 2]
  c23 <- costs[2, 3]
  c31 <- costs[3, 1]
  c32 <- costs[3, 2]
  c33 <- costs[3, 3] 
  rho1 <- rho[1]
  rho2 <- rho[2]
  rho3 <- rho[3]
  ### derivative means T1
  d111 <- th1st3(k1, rho1, n1, n2, n3, c11, c12, Thr[1])
  d121 <- th1st3(k2, rho2, n1, n2, n3, c21, c22, Thr[1])
  d131 <- th1st3(k3, rho3, n1, n2, n3, c31, c32, Thr[1])
  ### derivative standard deviations T1
  d112 <- th2st3(k1, rho1, n1, n2, n3, c11, c12, Thr[1])
  d122 <- th2st3(k2, rho2, n1, n2, n3, c21, c22, Thr[1])
  d132 <- th2st3(k3, rho3, n1, n2, n3, c31, c32, Thr[1])
  ### derivative means T2
  d211 <- th1st3(k1, rho1, n1, n2, n3, c12, c13, Thr[2])
  d221 <- th1st3(k2, rho2, n1, n2, n3, c22, c23, Thr[2])
  d231 <- th1st3(k3, rho3, n1, n2, n3, c32, c33, Thr[2])
  ### derivative standard deviations T2
  d212 <- th2st3(k1, rho1, n1, n2, n3, c12, c13, Thr[2])
  d222 <- th2st3(k2, rho2, n1, n2, n3, c22, c23, Thr[2])
  d232 <- th2st3(k3, rho3, n1, n2, n3, c32, c33, Thr[2])

  par1.1 <- mean(k1)
  par2.1 <- mean(k2)
  par3.1 <- mean(k3)
  par1.2 <- sd(k1)
  par2.2 <- sd(k2)
  par3.2 <- sd(k3)
  V1 <- d111^2*(par1.2^2/n1)+d121^2*(par2.2^2/n2) + d131^2*(par3.2^2/n3) + d112^2*(par1.2^2/(2*n1)) + d122^2*(par2.2^2/(2*n2)) + d132^2*(par3.2^2/(2*n3))
  V2 <- d211^2*(par1.2^2/n1)+d221^2*(par2.2^2/n2) + d231^2*(par3.2^2/n3) + d212^2*(par1.2^2/(2*n1)) + d222^2*(par2.2^2/(2*n2)) + d232^2*(par3.2^2/(2*n3))
  return(c(V1,V2))
}

##########################################################################
########### Partial derivative of cost function with respect to mean
##########################################################################
th1st3 <- function(k1, rho1, n1, n2, n3, c11, c12, Th){
  par1.1 <- mean(k1)
  par1.2 <- sd(k1)
  der <- (n1 + n2 + n3)/sqrt(2*pi)*rho1*(c11-c12)*((Th-par1.1)/par1.2^3)*exp(-(Th-par1.1)^2/(2*par1.2^2))
  return(der)
}
##########################################################################
########### Partial derivative of cost function with respect to standard deviation
##########################################################################
th2st3 <- function(k1, rho1, n1, n2, n3, c11, c12, Th){
  par1.1 <- mean(k1)
  par1.2 <- sd(k1)
  der <- -(n1 + n2 + n3)/(sqrt(2*pi)*par1.2^2)*rho1*(c11-c12)*exp(-(Th-par1.1)^2/(2*par1.2^2))*(-1+((Th-par1.1)^2/par1.2^2))
  return(der)
}

##########################################################################
########### Variance estimator
##################################################################################
VarThr3 <- function(k1, k2, k3, rho, Thr, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE)){
  c11 <- costs[1, 1]
  c12 <- costs[1, 2]
  c13 <- costs[1, 3]
  c21 <- costs[2, 1]
  c22 <- costs[2, 2]
  c23 <- costs[2, 3]
  c31 <- costs[3, 1]
  c32 <- costs[3, 2]
  c33 <- costs[3, 3] 
  rho1 <- rho[1]
  rho2 <- rho[2]
  rho3 <- rho[3]
  vd <- VDel3(k1, k2, k3, rho, Thr, costs)
  vn <- secondDer3aux(k1, k2, k3, rho, Thr, costs)
  ret1 <- vd[1]/vn[1]^2
  ret2 <- vd[2]/vn[2]^2
  # result
  re <- list(VAR1=ret1,VAR2=ret2)
  return(re)
}

#########################################################################
#### Derivatives of the cost function
#########################################################################
Der1 <- function(p, dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, rho, costs) {
  rho1 <- rho[1]
  rho2 <- rho[2]
  rho3 <- rho[3]
  c11 <- costs[1,1]
  c12 <- costs[1,2]
  c13 <- costs[1,3]
  c21 <- costs[2,1]
  c22 <- costs[2,2]
  c23 <- costs[2,3]
  c31 <- costs[3,1]
  c32 <- costs[3,2]
  c33 <- costs[3,3] 
  eq <- (rho1*(c11-c12)*dens(dist1)(p, par1.1, par1.2)+rho2*(c21-c22)*dens(dist2)(p, par2.1, par2.2)+rho3*(c31-c32)*dens(dist3)(p, par3.1, par3.2))
  return(eq)
}

Der2 <- function(p, dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, rho, costs) {
  rho1 <- rho[1]
  rho2 <- rho[2]
  rho3 <- rho[3]
  c11 <- costs[1, 1]
  c12 <- costs[1, 2]
  c13 <- costs[1, 3]
  c21 <- costs[2, 1]
  c22 <- costs[2, 2]
  c23 <- costs[2, 3]
  c31 <- costs[3, 1]
  c32 <- costs[3, 2]
  c33 <- costs[3, 3] 
  eq <- (rho1*(c12-c13)*dens(dist1)(p, par1.1, par1.2)+rho2*(c22-c23)*dens(dist2)(p, par2.1, par2.2)+rho3*(c32-c33)*dens(dist3)(p, par3.1, par3.2))
  return(eq)
}

##########################################################################
##### "uniroot" function looks for the one-variable equation solution  
##### arguments: 	q11, q12=percentiles of the first distribution 
#####			        q31, q32=percentiles of the third distribution 
#####		         dist1, dist2, dist3=choose the distributions
#####             between the following 2-parameter distributions: "beta", "cauchy",
#####             "chisq" (chi-squared), "gamma", "lnorm" (lognormal), "logis" (logistic), "norm" (normal)
#####             and "weibull".
#####			        par1.1=first parameter of the first distribution
#####		   	      par1.2=second parameter of the first distribution
#####             par2.1=first parameter of the second distribution
#####		          par2.2=second parameter of the second distribution
#####			        par3.1=first parameter of the third distribution
#####		          par3.2=second parameter of the third distribution
#####		          rho=vector of prevalences
#####		  	      costs=cost matrix
#####             tol=tolerance to be used in uniroot. Default, 10^(-8).
##### value: an object of class "thresTH3" containing the two thresholds,
#####        the prevalences and the cost matrix.  		  	
##################################################################################
thresTH3 <- function(dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE), q1=0.05, q2=0.5, q3=0.95, tol=10^(-8)){
  # error handling
  if (sum(rho > 0 & rho < 1) != 3){
    stop("The prevalences must be in (0,1)")
  }
  if (sum(rho) != 1){
    stop("The sum of the prevalences must be 1")
  }
  if (!is.matrix(costs)){
    stop("'costs' must be a matrix")
  }
  if (dim(costs)[1] != 3 | dim(costs)[2] != 3){
    stop("'costs' must be a 3x3 matrix")
  }
  costs.origin <- costs
  rho.origin <- rho
  # sort distributions
  median1 <- quant(dist1)(0.5, par1.1, par1.2)
  median2 <- quant(dist2)(0.5, par2.1, par2.2)
  median3 <- quant(dist3)(0.5, par3.1, par3.2)  
  median <- c(median1, median2, median3)
  ord_m <- order(median)
  if (!all(ord_m==1:3)){
    pars <- matrix(c(par1.1, par1.2, par2.1, par2.2, par3.1, par3.2), byrow=T, ncol=2)
    pars <- pars[ord_m, ]
    rho <- rho[ord_m]
    costs <- costs[ord_m, ord_m]
    par1.1 <- pars[1, 1]
    par1.2 <- pars[1, 2]
    par2.1 <- pars[2, 1]
    par2.2 <- pars[2, 2]
    par3.1 <- pars[3, 1]
    par3.2 <- pars[3, 2]
    dist <- c(dist1, dist2, dist3)[ord_m]
    dist1 <- dist[1]
    dist2 <- dist[2]
    dist3 <- dist[3]
  }
  # end of sorting
  p1 <- quant(dist1)(q1, par1.1, par1.2)
  p2 <- quant(dist2)(q2, par2.1, par2.2)
  p3 <- quant(dist3)(q3, par3.1, par3.2)
  cut1 <- uniroot(Der1, c(p1, p2), tol=tol, dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, rho, costs)$root
  cut2 <- uniroot(Der2, c(p2, p3), tol=tol, dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, rho, costs)$root
  # results
  re <- list(thres1=cut1, thres2=cut2, prev=rho.origin, costs=costs.origin, method="theoretical")
  class(re) <- "thresTH3"
  return(re)
}

############################################################################
# Print function for class "thresTH3"
############################################################################
print.thresTH3 <- function(x, ...){
    cat("\nThreshold 1: ", x$thres1)
    cat("\nThreshold 2: ", x$thres2)
    cat("\n")
    cat("\nParameters used")
    cat("\n  Prevalences:", x$prev)
    cat("\n  Costs")
    cat("\n    C11,C12,C13:", x$costs[1,])
    cat("\n    C21,C22,C23:", x$costs[2,])
    cat("\n    C31,C32,C33:", x$costs[3,])
    cat("\n")
}

##########################################################################
### cost function
##################################################################################
cost3fun <- function(x, k1, k2, k3, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE)){
  rho1 <- rho[1]
  rho2 <- rho[2]
  rho3 <- rho[3]
  c11 <- costs[1, 1]
  c12 <- costs[1, 2]
  c13 <- costs[1, 3]
  c21 <- costs[2, 1]
  c22 <- costs[2, 2]
  c23 <- costs[2, 3]
  c31 <- costs[3, 1]
  c32 <- costs[3, 2]
  c33 <- costs[3, 3] 
  n1 <- length(k1)
  n2 <- length(k2)
  n3 <- length(k3)
  ret <- (n1 + n2 + n3)*(rho1*(c11*pnorm(x[1], mean(k1), sd(k1)) + c12*(pnorm(x[2], mean(k1), sd(k1)) - pnorm(x[1], mean(k1), sd(k1))) + c13*(1-pnorm(x[2], mean(k1), sd(k1)))) + rho2*(c21*pnorm(x[1], mean(k2), sd(k2)) + c22*(pnorm(x[2], mean(k2), sd(k2)) - pnorm(x[1], mean(k2), sd(k2))) + c23*(1-pnorm(x[2], mean(k2), sd(k2)))) + rho3*(c31*pnorm(x[1], mean(k3), sd(k3)) + c32*(pnorm(x[2], mean(k3), sd(k3)) - pnorm(x[1], mean(k3), sd(k3))) + c33*(1-pnorm(x[2], mean(k3), sd(k3)))))
  return(ret)
}


################################################################################
############      DATA GENERATION
############ 	NON-PARAMETRIC RESAMPLING
##### arguments:   k1=vector containing the first sample values 
#####		   	       k2=vector containing the second sample values 
#####		   	       k3=vector containing the third sample values 
#####		   	       B=number of bootstrap resamples
##### returns: two-object list [[1]]:first resample matrix, [[2]]:second resample matrix, [[3]]:third resample matrix
################################################################################
resample3 <- function(k1, k2, k3, B){
  n1 <- length(k1)
  n2 <- length(k2)
  n3 <- length(k3)
  t1 <- matrix(sample(k1, n1*B, replace=TRUE), nrow=n1)
  t2 <- matrix(sample(k2, n2*B, replace=TRUE), nrow=n2)
  t3 <- matrix(sample(k3, n3*B, replace=TRUE), nrow=n3)
  t <- list(t1, t2, t3)
  return(t)
}


########################################################################
### 			Threshold estimator
#####		  	costs=cost matrix
########################################################################
thresNLM3 <- function(start, k1, k2, k3, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2],0 ,rho[3]/rho[2], 1, 1, 0), 3, 3,byrow=TRUE)){
  costs.origin <- costs
  rho.origin <- rho
  k1.origin <- k1
  k2.origin <- k2
  k3.origin <- k3
  # sort distributions 
  mean <- c(mean(k1), mean(k2), mean(k3))
  ord_m <- order(mean)
  if (!all(ord_m==1:3)){
    list.ks.ord <- list(k1=k1, k2=k2, k3=k3)[ord_m]
    k1 <- list.ks.ord[[1]]
    k2 <- list.ks.ord[[2]]
    k3 <- list.ks.ord[[3]]
    rho <- rho[ord_m]
    costs <- costs[ord_m, ord_m]
  }
  # end of sorting
  # estimate
  th <- nlm(cost3fun, start, k1, k2, k3, rho, costs)$estimate
  # results
  re <- list(thres1=th[1], thres2=th[2], prev=rho.origin, costs=costs.origin, k1=k1.origin, k2=k2.origin, k3=k3.origin)  
  return(re)
}

########################################################################
### 			Parametric Confidence Intervals
########################################################################
icParam3 <- function(k1, k2, k3, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0),3,3,byrow=TRUE), Thres, a=0.05){
  cut <- c(Thres[1],Thres[2])
  se0 <- VarThr3(k1,k2,k3,rho,cut,costs)
  se <- sqrt(c(se0$VAR1,se0$VAR2))
  # CIs
  ic1 <- c(cut[1]+qnorm(a/2)*se[1],cut[1]+qnorm(1-a/2)*se[1])
  ic2 <- c(cut[2]+qnorm(a/2)*se[2],cut[2]+qnorm(1-a/2)*se[2])
  # results
  ic <- list(lower1=ic1[1],upper1=ic1[2],lower2=ic2[1],upper2=ic2[2],alpha=a, ci.method="param")
  return(ic)
}



##########################################################################
#########         THRESHOLD COMPUTATION
##########################################################################
##### arguments: k1=vector containing the first sample values 
#####  	   	     k2=vector containing the second sample values
#####  	   	     k3=vector containing the third sample values
#####  	   	     rho=3-dimensional vector of prevalences
#####		   	     costs=cost matrix
#####            dist1, dist2, dist3: distribution to be assumed for the three populations,
#####            respectively. They can be chosen between the following 2-parameter distributions:
#####              "beta", "cauchy", "chisq" (chi-squared), "gamma", "lnorm" (lognormal),
#####              "logis" (logistic), "norm" (normal) and "weibull".
#####            ci.method=method to be used for the confidence intervals computation. The user can
#####            choose between:
#####              "param": parametric CIs are computed when assuming a trinormal underlying model.
#####              "boot": the confidence interval is computed by bootstrap.
#####            Default, "param". The user can specify just just the initial letters.
#####            B=number of bootstrap resamples when ci.method="boot". Otherwise, ignored.
#####            Default, 1000.
#####  	   	     alpha=significance level for the confidence interval. Default, 0.05.
#####            na.rm=a logical value indicating whether NA values in k1, k2 and k3 should
#####            be stripped before the computation proceeds. Default, F.
##### value: the thresholds estimated
#######################################
thres3 <- function(k1, k2, k3, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2],0,rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE), dist1="norm", dist2="norm", dist3="norm", start=NULL, ci.method=c("param", "boot"), B=1000, alpha=0.05, na.rm=FALSE){
  # error handling
  if (sum(rho > 0 & rho < 1) != 3){
    stop("The prevalences must be in (0,1)")
  }
  if (sum(rho) != 1){
    stop("The sum of the prevalences must be 1")
  }
  if (!is.matrix(costs)){
    stop("'costs' must be a matrix")
  }
  if (dim(costs)[1] != 3 | dim(costs)[2] != 3){
    stop("'costs' must be a 3x3 matrix")
  }
  if (!is.numeric(k1) | !is.numeric(k2) | !is.numeric(k3)){
    stop("'k1', 'k2' and 'k3' must be numeric vectors")
  }
  # NAs handling
  if (na.rm){
    k1 <- k1[!is.na(k1)]
    k2 <- k2[!is.na(k2)]
    k3 <- k3[!is.na(k3)]
  }
  # function
  ci.method <- match.arg(ci.method)
  if (dist1=="norm" & dist2=="norm" & dist3=="norm"){
    if (is.null(start)){
      stop("'start' must be specified")
    }
    if (length(start)!=2){
      stop("'start' must be a 2-dimensional vector containing starting values for the thresholds") 
    }
    T <- thresNLM3(start, k1, k2, k3, rho, costs)
    T$dist1 <- dist1
    T$dist2 <- dist2
    T$dist3 <- dist3
    if (ci.method=="param"){
      ci <- icParam3(k1, k2, k3, rho, costs, c(T$thres1, T$thres2), a=alpha)
    }
    if (ci.method=="boot"){
      ci <- icBoot3(k1, k2, k3, rho, costs, c(T$thres1, T$thres2), B=B, a=alpha, start)
    }
  }else{
    if (ci.method=="param"){
      stop("When not all the distributions are 'norm', parametric CIs cannot be computed (choose ci.method='boot')")
    }
    if (!(dist1 %in% c("beta", "cauchy", "chisq", "gamma", "lnorm", "logis", "nbinom", "norm", "weibull"))){
      stop("Unsupported distribution for 'dist1'")
    }
    if (!(dist2 %in% c("beta", "cauchy", "chisq", "gamma", "lnorm", "logis", "nbinom", "norm", "weibull"))){
      stop("Unsupported distribution for 'dist2'")
    }
    if (!(dist3 %in% c("beta", "cauchy", "chisq", "gamma", "lnorm", "logis", "nbinom", "norm", "weibull"))){
      stop("Unsupported distribution for 'dist3'")
    }
    # parameter estimation through 'fitdistr()'
    # dist1
    pars1 <- getParams(k1, dist1)
    # dist2
    pars2 <- getParams(k2, dist2)
    # dist3
    pars3 <- getParams(k3, dist3)
    # threshold+CI estimation
    T <- thresTH3(dist1, dist2, dist3, pars1[1], pars1[2], pars2[1], pars2[2], pars3[1], pars3[2], rho, costs, tol=10^(-8))
    T <- unclass(T)
    # adding some information to the output
    T$k1 <- k1
    T$k2 <- k2
    T$k3 <- k3
    T$dist1 <- dist1
    T$dist2 <- dist2
    T$dist3 <- dist3
    T$pars1 <- pars1
    T$pars2 <- pars2
    T$pars3 <- pars3
    # confidence interval by parametric bootstrap
    ci <- icBootTH3(dist1, dist2, dist3, pars1[1], pars1[2], pars2[1], pars2[2], pars3[1], pars3[2], length(k1), length(k2), length(k3), rho, costs, c(T$thres1, T$thres2), B, alpha)
  }  
  out <- list(T=T, CI=ci)
  class(out) <- "thres3"
  return(out)
}



############################################################################
# Print function for class "thres3"
############################################################################
print.thres3 <- function(x, ...){
  if (x$T$dist1 == "norm" & x$T$dist2 == "norm" & x$T$dist3 == "norm"){
    cat("\nEstimate:")
    cat("\n  Threshold 1: ", x$T$thres1)
    cat("\n  Threshold 2: ", x$T$thres2)
    cat("\n")
    if(x$CI$ci.method == "param"){
      cat("\nConfidence interval Threshold 1:")
      cat("\n  Lower Limit:", x$CI$lower1)
      cat("\n  Upper Limit:", x$CI$upper1)
      cat("\n")
      cat("\nConfidence interval Threshold 2:")
      cat("\n  Lower Limit:", x$CI$lower2)
      cat("\n  Upper Limit:", x$CI$upper2)
      cat("\n")
    }
    if(x$CI$ci.method == "boot"){
      cat("\nConfidence intervals (bootstrap):")
      cat("\n  CI based on normal distribution for Threshold 1: ", x$CI$low.norm1, " - ", x$CI$up.norm1)
      cat("\n  CI based on percentiles for Threshold 1: ", x$CI$low.perc1, " - ", x$CI$up.perc1)
      cat("\n  CI based on normal distribution for Threshold 2: ", x$CI$low.norm2, " - ", x$CI$up.norm2)
      cat("\n  CI based on percentiles for Threshold 2: ", x$CI$low.perc2, " - ", x$CI$up.perc2)
      cat("\n  Bootstrap resamples: ", x$CI$B)
      cat("\n")
    }    
    cat("\nParameters used:")
    cat("\n  Prevalences:", x$T$prev)
    cat("\n  Costs")
    cat("\n    C11,C12,C13:", x$T$costs[1,])
    cat("\n    C21,C22,C23:", x$T$costs[2,])
    cat("\n    C31,C32,C33:", x$T$costs[3,])
    cat("\n  Significance Level: ", x$CI$alpha)
    cat("\n")

  }else{
    cat("\nEstimate:")
    cat("\n  Threshold 1: ", x$T$thres1)
    cat("\n  Threshold 2: ", x$T$thres2)
    cat("\n")
    cat("\nConfidence intervals (parametric bootstrap):")
      cat("\n  CI based on normal distribution for Threshold 1: ", x$CI$low.norm1, " - ", x$CI$up.norm1)
      cat("\n  CI based on percentiles for Threshold 1: ", x$CI$low.perc1, " - ", x$CI$up.perc1)
      cat("\n  CI based on normal distribution for Threshold 2: ", x$CI$low.norm2, " - ", x$CI$up.norm2)
      cat("\n  CI based on percentiles for Threshold 2: ", x$CI$low.perc2, " - ", x$CI$up.perc2)
      cat("\n  Bootstrap resamples: ", x$CI$B)
      cat("\n")
    cat("\nParameters used:")
    cat("\n  Prevalences:", x$T$prev)
    cat("\n  Costs")
    cat("\n    C11,C12,C13:", x$T$costs[1,])
    cat("\n    C21,C22,C23:", x$T$costs[2,])
    cat("\n    C31,C32,C33:", x$T$costs[3,])
    cat("\n  Confidence Level: ", x$CI$alpha)
    cat("\n  Distribution assumed for the first sample: ", x$T$dist1, "(", round(x$T$pars1[1], 2), ", ", round(x$T$pars1[2], 2), ")", sep="")
    cat("\n  Distribution assumed for the second sample: ", x$T$dist2, "(", round(x$T$pars2[1], 2), ", ", round(x$T$pars2[2], 2), ")", sep="")
    cat("\n  Distribution assumed for the third sample: ", x$T$dist3, "(", round(x$T$pars3[1], 2), ", ", round(x$T$pars3[2], 2), ")", sep="") 
    cat("\n")
  }
}

########################################################################
### 			Bootstrap Confidence Interval
########################################################################
icBoot3 <- function(k1, k2, k3, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE), Thres, B, a=0.05, start){
  t <- resample3(k1, k2, k3, B)
  t1 <- t[[1]]
  t2 <- t[[2]]
  t3 <- t[[3]]
  # bootstrap
  cutsim <- sapply(1:B,function(i){
    thresNLM3(start,t1[,i],t2[,i],t3[,i],rho,costs)[1:2]
  })
  
  cutsim1 <- unlist(cutsim[1, ])
  cutsim2 <- unlist(cutsim[2, ])
  
  # sd
  est.se1 <- sd(cutsim1)
  est.se2 <- sd(cutsim2)
  
  
  ###### 1) NORMAL-BOOTSTRAP SE
  norm1 <- c(Thres[1] + qnorm(a/2)*est.se1, Thres[1] + qnorm(1-a/2)*est.se1)
  norm2 <- c(Thres[2] + qnorm(a/2)*est.se2, Thres[2] + qnorm(1-a/2)*est.se2)
  
  ###### 2) PERCENTILE
  perc1 <- c(quantile(cutsim1, a/2), quantile(cutsim1, 1-a/2))
  perc2 <- c(quantile(cutsim2, a/2), quantile(cutsim2, 1-a/2))
  
  # results
  re <- list(low.norm1=norm1[1], up.norm1=norm1[2], low.norm2=norm2[1], up.norm2=norm2[2], low.perc1=perc1[1], up.perc1=perc1[2], low.perc2=perc2[1], up.perc2=perc2[2], alpha=a, B=B, ci.method="boot")
  return(re)
}

################################################################################
############       DATA GENERATION
############    PARAMETRIC RESAMPLING
##### arguments:  dist1, dist2, dist3=distributions assumed for the three populations. 
#####  	   	      B=number of bootstrap resamples
##### value: three-object list [[1]]:first population resample matrix, [[2]]:second population resample matrix,
#####         [[3]]: third population resample matrix
################################################################################
aux.par.boot3 <- function(dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, n1, n2, n3, B){
  t0 <- matrix(rand(dist1)(n1*B, par1.1, par1.2), nrow=n1) # B resamples of dist1(par1.1, par1.2) with replacement
  t1 <- matrix(rand(dist2)(n2*B, par2.1, par2.2), nrow=n2) # B resamples of dist2(par2.1, par2.2) with replacement
  t2 <- matrix(rand(dist3)(n3*B, par3.1, par3.2), nrow=n3) # B resamples of dist3(par3.1, par3.2) with replacement
  t <- list(t0, t1, t2)
  return(t)
}


##########################################################################
#########  GETTING PARAMETERS OF A DISTRIBUTION THROUGH fitdistr() [library MASS]
##########################################################################
getParams <- function(k, dist){
  if (dist %in% c("cauchy", "gamma", "weibull")){
    pars <- fitdistr(k, dist)$estimate
  }else if(dist=="beta"){
    # needs initial values
    sigma2 <- var(k)
    mu <- mean(k)
    shape1.start <- ((1-mu)/sigma2-1/mu)*mu^2
    shape2.start <- shape1.start*(1/mu-1)
    pars <- fitdistr(k, "beta", start=list(shape1=shape1.start, shape2=shape2.start))$estimate
  }else if(dist=="chisq"){
    # needs initial values
    sigma2 <- var(k)
    mu <- mean(k)
    ncp.start <- sigma2/2-mu
    df.start <- mu-ncp.start
    pars <- fitdistr(k, "chi-squared", start=list(df=df.start, ncp=ncp.start))$estimate
  }else if(dist=="lnorm"){
    pars <- fitdistr(k, "lognormal")$estimate
  }else if(dist=="logis"){
    pars <- fitdistr(k, "logistic")$estimate
  }else if(dist=="norm"){
    pars <- fitdistr(k, "normal")$estimate
  }
  out <- pars
  return(out)
}


###################################################################################
###################################################################################
###############    FUNCTION OF PARAMETRIC BOOTSTRAP FOR THEORETICAL METHOD
##########              1) BOOTSTRAP VERSION OF S.E.- NORMAL
##########              2) PERCENTILE BIAS-CORRECTED
##### arguments:  dist1, dist2, dist3 =distributions assumed for three populations. 
#####             par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, n1, n2, n3 =parameters and sample sizes
#####		   	      B=number of bootstrap resamples
#####		   	      a=significance level
#####		   	      rho=disease prevalences
#####		  	      costs=cost matrix
##### returns: the 2 interval limits
##################################################################################
icBootTH3 <- function(dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, n1, n2, n3, rho, costs=matrix(c(0, 1, 1, rho[1]/rho[2], 0, rho[3]/rho[2], 1, 1, 0), 3, 3, byrow=TRUE), Thres, B, a){
  # bootstrap resamples
  t <- aux.par.boot3(dist1, dist2, dist3, par1.1, par1.2, par2.1, par2.2, par3.1, par3.2, n1, n2, n3, B)
  t0 <- t[[1]]
  t1 <- t[[2]]
  t2 <- t[[3]]
  # computing threshold of the resamples...
  # 1) parameter estimation for each resample
  pars1 <- sapply(1:B, function(i){getParams(t0[, i], dist1)})
  pars2 <- sapply(1:B, function(i){getParams(t1[, i], dist2)})
  pars3 <- sapply(1:B, function(i){getParams(t2[, i], dist3)})
  # 2) TH cut
  cut <- sapply(1:B, function(i){thresTH3(dist1, dist2, dist3, pars1[1, i], pars1[2, i], pars2[1, i], pars2[2, i], pars3[1,i], pars3[2,i], rho, costs)[1:2]})
  
  cut1 <- unlist(cut[1, ])
  cut2 <- unlist(cut[2, ])

  est.se1 <- sd(cut1)
  est.se2 <- sd(cut2)
  
  ###### 1) NORMAL-BOOTSTRAP SE
  norm1 <- c(Thres[1] + qnorm(a/2)*est.se1, Thres[1] + qnorm(1-a/2)*est.se1)
  norm2 <- c(Thres[2] + qnorm(a/2)*est.se2, Thres[2] + qnorm(1-a/2)*est.se2)
  
  ###### 2) PERCENTILE
  perc1 <- c(quantile(cut1, a/2), quantile(cut1, 1-a/2))
  perc2 <- c(quantile(cut2, a/2), quantile(cut2, 1-a/2))
  
  # results
  re <- list(low.norm1=norm1[1], up.norm1=norm1[2], low.norm2=norm2[1], up.norm2=norm2[2], low.perc1=perc1[1], up.perc1=perc1[2], low.perc2=perc2[1], up.perc2=perc2[2], alpha=a, B=B, ci.method="boot")
  return(re)
}



##################################################################################
######	DENSITY PLOT
##################################################################################
###### plot.thres3 function provides a graphic including the three sample densities,
###### the thresholds and their confidence intervals
###### arguments:
######  x: 'thres3' object
######  bw: vector containing the bandwith for the 1st sample in the 1st
######      position, the bandwith for the 2nd sample in the 2nd position
######      and bandwith for the 3rd sample in the 3rd position
######      (to be passed to 'density()'). Default, c("nrd0", "nrd0", "nrd0").
######  ci: should the confidence intervals be plotted? Default, T.
######  which.boot: in case 'x' contains CI computed by bootstrapping, which one should be printed?
######              the user can choose between "norm" (based on normal distribution)
######              or "perc" (based on percentiles). Default, "norm". This argument
######              is ignored if parametric CI were computed 
######  col: 4-dimensional vector:
######       col[1]: color for the density of the first sample
######       col[2]: color for the density for the second sample
######       col[3]: color for the density of the third sample
######       col[4]: color for the thresholds and their corresponding CI
######       Default, c(1, 2, 3, 1).
######  lty: 5-dimensional vector:
######       lty[1]: line type for the density of the first sample
######       lty[2]: line type for the density of the second sample
######       lty[3]: line type for the density of the third sample
######       lty[4]: line type for the thresholds
######       lty[5]: line type for the CI
######       Default, c(1, 1, 1, 1, 2).
######  lwd: 4-dimensional vector:
######       lwd[1]: line width for the density of the first sample
######       lwd[2]: line width for the density for the second sample
######       lwd[3]: line width for the density for the third sample
######       lwd[4]: line width for the thresholds and their corresponding CI
######       Default, c(1, 1, 1, 1).
######  main, xlab, ...: further arguments to be passed to 'plot()'.
######  legend: logical asking if an automatic legend should be added to the graph. Default, TRUE.
######  leg.pos: position of the legend. Default, "topleft". Ignored if legend=FALSE.
######  leg.cex: a number that reescales the size of the legend. Ignored if legend=FALSE. Default, 1.
##################################################################################
plot.thres3 <- function(x, bw=c("nrd0", "nrd0", "nrd0"), ci=TRUE, which.boot=c("norm", "perc"), col=c(1, 2, 3, 1), lty=c(1, 1, 1, 1, 2), lwd=c(1, 1, 1, 1), main=paste0("Threshold estimates", ifelse(ci, " and CIs", "")), xlab="", legend=TRUE, leg.pos="topleft", leg.cex=1, ...){
  which.boot <- match.arg(which.boot)
  k1 <- x$T$k1
  k2 <- x$T$k2
  k3 <- x$T$k3
  dens.k1 <- density(k1, bw=bw[1])
  dens.k2 <- density(k2, bw=bw[2])
  dens.k3 <- density(k3, bw=bw[3])
  min.x <- min(min(dens.k1$x), min(dens.k2$x), min(dens.k3$x))
  max.x <- max(max(dens.k1$x), max(dens.k2$x), max(dens.k3$x))
  max.y <- max(max(dens.k1$y), max(dens.k2$y), max(dens.k3$y))
  plot(dens.k1, xlim=c(min.x, max.x), ylim=c(0, max.y), col=col[1], lty=lty[1], lwd=lwd[1], main=main, xlab=xlab, ...)
  lines(dens.k2, col=col[2], lty=lty[2], lwd=lwd[2])
  lines(dens.k3, col=col[3], lty=lty[3], lwd=lwd[3])
  # thres
  abline(v=x$T$thres1, col=col[4], lty=lty[4], lwd=lwd[4])
  abline(v=x$T$thres2, col=col[4], lty=lty[4], lwd=lwd[4])
  # CI
  if(ci){
    if (x$CI$ci.method != "boot"){
      abline(v=c(x$CI$lower1, x$CI$upper1), col=col[4], lty=lty[5], lwd=lwd[4])
      abline(v=c(x$CI$lower2, x$CI$upper2), col=col[4], lty=lty[5], lwd=lwd[4])
    }else{
      abline(v=c(x$CI[paste0("low.", which.boot, "1")], x$CI[paste0("up.", which.boot, "1")]), col=col[4],  lty=lty[5], lwd=lwd[4])
      abline(v=c(x$CI[paste0("low.", which.boot, "2")], x$CI[paste0("up.", which.boot, "2")]), col=col[4],  lty=lty[5], lwd=lwd[4])
    }
  }  
  # legend
  if (legend){
    legend(leg.pos, c("1st sample", "2nd sample", "3rd sample", ifelse(ci, "Thres+CI", "Thres")), col=col, lty=lty, lwd=lwd, cex=leg.cex)
  }
}



##################################################################################
######  LINES DENSITY PLOT
##################################################################################
###### lines.thres3 function includes vertical lines for the thresholds and CIs in
###### a plot.thres3.
###### arguments:
######  x: 'thres3' object
######  ci: should the confidence intervals be plotted? Default, TRUE.
######  which.boot: in case 'x' contains CI computed by bootstrapping, which one should be printed?
######              the user can choose between "norm" (based on normal distribution)
######              or "perc" (based on percentiles). Default, "norm". This argument
######              is ignored if parametric CI were computed 
######  col: color for the thresholds and their corresponding CI. Default, 1.
######  lty: 2-dimensional vector:
######       lty[1]: line type for the thresholds
######       lty[2]: line type for the CI
######       Default, c(1, 2).
######  lwd: line width for the thresholds and their corresponding CI. Default, 1.
######  ...: further arguments to be passed to 'abline()'.
##################################################################################
lines.thres3 <- function(x, ci=TRUE, which.boot=c("norm", "perc"), col=1, lty=c(1, 2), lwd=1, ...){
  which.boot <- match.arg(which.boot)
  # thres
  abline(v=x$T$thres1, col=col, lty=lty[1], lwd=lwd)
  abline(v=x$T$thres2, col=col, lty=lty[1], lwd=lwd)
  # CI
  if(ci){
    if (x$CI$ci.method != "boot"){
      abline(v=c(x$CI$lower1, x$CI$upper1), col=col, lty=lty[2], lwd=lwd)
      abline(v=c(x$CI$lower2, x$CI$upper2), col=col, lty=lty[2], lwd=lwd)
    }else{
      abline(v=c(x$CI[paste0("low.", which.boot, "1")], x$CI[paste0("up.", which.boot, "1")]), col=col,  lty=lty[2], lwd=lwd)
      abline(v=c(x$CI[paste0("low.", which.boot, "2")], x$CI[paste0("up.", which.boot, "2")]), col=col,  lty=lty[2], lwd=lwd)
    }
  }  
}
"DeltaMethod" <- function(par, fct, vcov, delta, ...){
#
# DeltaMethod  - computes delta method v-c matrix of results of function fct at values of par
#                using var-cov matrix (vcov) of parameters and numerical first derivatives 
#                of fct with respect to par.  
#
# Arguments:
#
# par - parameter estimates 
# fct - function that produces estimate
# vcov - var-cov matrix of par
# delta - relative change in par to compute central difference first derivative
# ...   - any remaining additional arguments needed for fct
#
# Value:
#
# variance-cov matrix of estimates 
#
#   Construct theta call to fct
#
  theta=function(par) fct(par,...)

  #   Numerically compute the first derivative of the estimator with respect to each parameter
  #   using a central difference formula
  savepar<-par
  value1=theta(par)
  partial<-matrix(0,nrow=length(par),ncol=length(value1))
  for(i in 1:length(par)){
    # Store the original parameters into par and then adjust the ith parameter by adding the proportion 
    # based on delta
    par<-savepar
    par[i]<-savepar[i]+delta*savepar[i] 
    # With this new value call the function to compute the estimate
    value1=theta(par)
    # Next do the same thing as above except substract off the proportion delta from the original parameter.
    par<-savepar
    par[i]<-savepar[i]-delta*savepar[i]    
    value2=theta(par)
    # Compute the central difference formula for the first derivative with respect to the ith parameter
    partial[i,]<-(value1-value2)/(2*delta*savepar[i])
  }

  variance<-t(partial)%*%vcov%*%partial 
#
#  1/26/06 jll; modified to return the v-c matrix and the first partial vector(matrix)
#
  return(list(variance=variance,partial=partial))
}

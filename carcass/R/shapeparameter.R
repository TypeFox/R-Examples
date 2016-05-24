shapeparameter <-
function(m, lwr=NA, upr=NA, se=NA){
# m = estimate
# lwr, upr = lower and upper limit of the 95 % credible or confidence 
# se = standard error of the proportion
# interval
#--------------------------------------------------------------------
if(is.na(se)){
  if(sum(lwr<0)>0) stop("Limits of confidence interval need to be between 0 and 1.")
  if(sum(upr>1)>0) stop("Limits of confidence interval need to be between 0 and 1.")
  if(sum(upr<lwr)>0) stop("Lower limit of confidence interval needs to be lower than the upper limit.")
  ci <- upr - lwr
  sigma2 <-(ci/4)^2
  }
if(!is.na(se)){
  sigma2 <- se^2
  }
a <- m*(m*(1-m)/sigma2-1)
b <- (1-m)*(m*(1-m)/sigma2-1)
if(sum(a<0)>0|sum(b<0)>0) stop("A too large variance for a probability parameter has been specified.")
return(list(a=a,b=b))
}


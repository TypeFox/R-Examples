LRNormalMean_KV <-
function(x,mu0,alpha,sigma)  {
  ifelse(abs(sqrt(length(x))*(mean(x)-mu0)/sigma)>qnorm(1-alpha/2),"Reject Hypothesis H","Fail to Reject Hypothesis H")
}

LRNormalMean_UV <-
function(x,mu0,alpha){
  S <- sd(x); n <- length(x)
  ifelse(abs(sqrt(length(x))*(mean(x)-mu0)/S)>qt(n-1,1-alpha/2),"Reject Hypothesis H","Fail to Reject Hypothesis H")
}

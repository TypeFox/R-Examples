rgammapois <-
function(n, p, r){
  choosefrom <- 0:10000
  probs <- exp(lgamma(r+choosefrom)-lgamma(choosefrom+1)-lgamma(r)+choosefrom*log(p)+r*log(1-p))
  mysample <- sample(x=choosefrom, size=n, replace=T, prob=probs)
  return(mysample)
}

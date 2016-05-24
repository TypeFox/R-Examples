dpolya <-
function(x=NA,alpha=NA){
  nfac<-factorial(sum(x))
  pnfac<-prod(factorial(x))
  term1<-nfac / pnfac
  sumAlpha<-sum(alpha)
  term2<-factorial(sumAlpha - 1) / factorial(sum(x) + sumAlpha - 1)
  term3<-prod(factorial(x + alpha - 1) / factorial(alpha - 1))
  prob<-term1 * term2 * term3
  return(prob)
}


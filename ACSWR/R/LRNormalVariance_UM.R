LRNormalVariance_UM <-
function(x,sigma0,alpha){
  S <- var(x); n <- length(x)
  chidata <- ((n-1)*S)/(sigma0^2)
  ifelse((chidata<qchisq(df=n-1,p=alpha/2)|| (chidata>qchisq(df=n-1,p=1-alpha/2))),"Reject Hypothesis H","Fail to Reject Hypothesis H")
}

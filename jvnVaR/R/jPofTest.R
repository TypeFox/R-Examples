jPofTest <-
function(n,k,p,test_significant){
# 
statistic <- -2*log(((1-p)^(n-k)*p^k)/((1-k/n)^(n-k)*(k/n)^k))
Quantile <- qchisq(1-test_significant,1)
rslt <- statistic <= Quantile
return(c(statistic,Quantile,rslt))
}

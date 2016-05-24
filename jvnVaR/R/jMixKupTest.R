jMixKupTest <-
function(n,v,p,test_significant){
# 
k <- length(v)
test <- jPofTest(n,k,p,test_significant)
statistic <- test[1]
for (i in 1:k) {
test <- jTuffTest(n,v[i],p,test_significant)
statistic <- statistic + test[1]
}
Quantile <- qchisq(1-test_significant,k+1)
rslt <- statistic <= Quantile
return(c(statistic,Quantile,rslt))
}

jTuffTest <-
function(n,v,p,test_significant){
# 
statistic <- -2*log(p*(1-p)^(v-1)/((1/v)*(1-1/v)^(v-1)));
Quantile <- qchisq(1-test_significant,1)
rslt <- statistic <= Quantile
return(c(statistic,Quantile,rslt))
}

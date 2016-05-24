testlik.fun <-
function(ModG,ModR)
{

difdev<-2*(-ModG@min+ModR@min)
gl<-ModG@npar-ModR@npar
pv<-1-pchisq(difdev,gl)
cat(fill=TRUE)
cat('General Model (hypothesis H1): ',ModG@tit, fill=T)
cat('Reduced Model (hypothesis H0): ',ModR@tit, fill=T)
cat('ML ratio test statistic: ', round(difdev,2),fill=TRUE)
cat('P-value:', round(pv,3),fill=TRUE)
cat(fill=TRUE)
return(list(pv=pv,ModG=ModG,ModR=ModR))
}

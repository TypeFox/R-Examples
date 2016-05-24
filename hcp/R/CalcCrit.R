CalcCrit <-
function(fit.select){sumObj<-summary(fit.select)
 R2ad<-sumObj$adj.r.squared
SSE<-tail(anova(fit.select)$"Sum Sq",1)
n<-length(fit.select$fitted)
p<-fit.select$rank
AIC<-n*log(SSE)-n*log(n)+2*p
BIC<-n*log(SSE)-n*log(n)+log(n)*p
return(c(R2ad,AIC,BIC))
}

corr.nn4bn<-function(p, BN.cor){
corrected=BN.cor/(dnorm(qnorm(p))/sqrt(p*(1-p)))
return(corrected)
}

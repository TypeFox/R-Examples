CItran.fun <-
function(VARbeta, lambdafit, covariates, clevel=0.95)
{
covariates<<-covariates
VARbeta<<-VARbeta
n<-dim(covariates)[1]
VarXbeta<-rep(NA,n)
calcVlambda.fun<-function(i)
{
VarXbeta[i]<<-matrix(covariates[i,],nrow=1)%*%VARbeta%*%matrix(covariates[i,],ncol=1)
}

aux<-sapply(c(1:n), FUN=calcVlambda.fun)

UIlambda<-exp(log(lambdafit)+ qnorm((1-(1-clevel)/2))*VarXbeta**0.5)
LIlambda<-exp(log(lambdafit)- qnorm((1-(1-clevel)/2))*VarXbeta**0.5)
return(list(UIlambda=UIlambda, LIlambda=LIlambda, lambdafit=lambdafit))
}



VARbeta.fun <-
function(covariates, lambdafit)
{

lambdafit<-lambdafit
covariates<-covariates

K<-dim(covariates)[2]
InforM<-matrix(rep(0,K*K), ncol=K)

calcVAR.fun<-function(j,i)
{
InforM[i,j]<<-sum(lambdafit*covariates[,i]*covariates[,j])
InforM[j,i]<<-sum(lambdafit*covariates[,i]*covariates[,j])
}

calcVAR2.fun<-function(i)
{
aux<-sapply(c(1:i), FUN=calcVAR.fun, i)
}

aux2<-sapply(c(1:K), FUN=calcVAR2.fun)
VARbeta<-try(solve(InforM))
if (is.matrix(VARbeta))  attr(VARbeta,'CalMethod')<-'solve function'
else 
{
VARbeta<-try(chol2inv(chol(InforM)))
if (is.matrix(VARbeta))  attr(VARbeta,'CalMethod')<-'Cholesky'
else {
	VARbeta<-matrix(numeric(), 0L, 0L)
	 attr(VARbeta,'CalMethod')<-'Not possible'
	}


}


return(VARbeta)

}

nntsmanifoldnewtonestimation <-
function(data,M=0,iter=1000,initialpoint=FALSE,cinitial)
{
	data<-as.matrix(data)
	n<-length(data)
	R<-1
	if (R != length(M))
		return("Error: M must correspond to a univariate variable distribution")

	statisticsmatrix<-matrix(0,nrow=M+1,ncol=n)

	for (k in 0:M)
		statisticsmatrix[k+1,]<-t(Conj(exp(1i*k*data)))

	if (initialpoint)
	{
		if (length(cinitial) != (M+1))
			return("Length of cinitial must be equal to M+1")
		c0<-cinitial
	}
	else
		c0<-apply(statisticsmatrix,1,mean)

	c0<-c0/sqrt(sum(Mod(c0)^2))

	eta<-matrix(0,nrow=M+1,ncol=1)

	for (k in 1:n)
		eta <- eta + (1/n)*(1/(t(Conj(c0))%*%statisticsmatrix[ ,k]))*statisticsmatrix[ ,k]

	eta<-eta-c0
	
	newtonmanifold<-(c0+eta)
	newtonmanifold<-newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
	newtonmanifold<-newtonmanifold*exp(-1i*Arg(newtonmanifold[1]))
	newtonmanifoldprevious<-newtonmanifold

	for (j in 1:iter){
		eta<-matrix(0,nrow=M+1,ncol=1)
		for (k in 1:n){
			eta<-eta+(1/n)*(1/(t(Conj(newtonmanifold))%*%statisticsmatrix[ ,k]))*statisticsmatrix[ ,k]
			}
		eta<-eta-newtonmanifold
		newtonmanifold<-newtonmanifold+eta
		newtonmanifold<-newtonmanifold/sqrt(sum(Mod(newtonmanifold)^2))
		newtonmanifold<-newtonmanifold*exp(-1i*Arg(newtonmanifold[1]))
		if (j==iter)
			normsequence<-(sqrt(sum(Mod(newtonmanifold-newtonmanifoldprevious)^2)))
		newtonmanifoldprevious<-newtonmanifold
		}

	newtonmanifold<-newtonmanifold/sqrt(2*pi)
	loglik <- nntsloglik(data,newtonmanifold,M) 
	AIC <- -2*loglik + 2*(2*M)
	BIC <- -2*loglik + (2*M)*log(n)
	gradnormerror <- normsequence

	cestimatesarray <- data.frame(cbind(0:M,newtonmanifold))
	cestimatesarray[,1]<-as.integer(Re(as.matrix(cestimatesarray[,1])))
	names(cestimatesarray)[1]<-"k"
	names(cestimatesarray)[2]<-"cestimates"

	res<-list(cestimates=cestimatesarray, loglik=loglik, AIC=AIC, BIC=BIC, gradnormerror=gradnormerror)
	return(res)
}


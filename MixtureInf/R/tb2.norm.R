tb2.norm <-
function(alpha,theta)
{
#This function computes \tilde B_{22} under normal mixture.
#\tilde B_{22} has been standardized to the correlation matrix.
#
#To calculate \tilde B_{22}, (-\infty,\infty) has been divided into N=10000 subintervals 
#with the probability in each subinterval being 1/N.
#
#In practice, you may change N to be a large value to see its effect.
#
#alpha:  vector of mixture probabilities.
#theta:  vector of means of each component.
	m0=length(alpha)
	N=10000
	quan=matrix((0:(N-1)+0.5)/N,ncol=1)
	x=as.numeric(apply(quan,1,qmixnorm,alpha=alpha,theta=theta))
	pdf=0
	delta=c()
	y=c()
	z=c()

	for(i in 1:m0)
	{
		delta=cbind(delta,dnorm(x,theta[i],1)-dnorm(x,theta[m0],1))
		y=cbind(y,(x-theta[i])*dnorm(x,theta[i],1))
		z=cbind(z,((x-theta[i])^2-1)*dnorm(x,theta[i],1)/2 )
		pdf=pdf+alpha[i]*dnorm(x,theta[i],1)
	}
	bi=cbind(delta[,1:(m0-1)],y,z)/pdf
	B=t(bi)%*%bi	
	B11=B[1:(2*m0-1),1:(2*m0-1)]
	B12=B[1:(2*m0-1),(2*m0):(3*m0-1)]
	B22=B[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]
	tB22=B22-t(B12)%*%solve(B11)%*%B12
	diagB22=diag(diag(tB22)^(-1/2),m0,m0)
	corr=diagB22%*%tB22%*%diagB22
	round(corr,3)
}

tb2.binom <-
function(alpha,theta,size)
#This function computes \tilde B_{22} under binomial mixture.
#
#\tilde B_{22} has been standardized to the correlation matrix.
#
#alpha:  vector of mixture probabilities.
#theta:  vector of probabilities of success of each component.
#size:   number of trials.
{
	m0=length(alpha)
	x=0:size
	pdf=0
	delta=c()
	y=c()
	z=c()

	for(i in 1:m0)
	{
		delta=cbind(delta,dbinom(x,size,theta[i])-dbinom(x,size,theta[m0]))
		d1binom=(x/theta[i]-(size-x)/(1-theta[i]))*dbinom(x,size,theta[i])
		y=cbind(y,d1binom)
		d2binom=(x*(x-1)/theta[i]^2-2*x*(size-x)/(theta[i]*(1-theta[i]))+(size-x)*(size-x-1)/((1-theta[i])^2))*dbinom(x,size,theta[i])
		z=cbind(z,d2binom/2)
		pdf=pdf+alpha[i]*dbinom(x,size,theta[i])
	}
	bi=cbind(delta[,1:(m0-1)],y,z)/pdf
	B=t(bi)%*%diag(pdf)%*%bi
	B11=B[1:(2*m0-1),1:(2*m0-1)]
	B12=B[1:(2*m0-1),(2*m0):(3*m0-1)]
	B22=B[(2*m0):(3*m0-1),(2*m0):(3*m0-1)]
	tB22=B22-t(B12)%*%solve(B11)%*%B12
	diagB22=diag(diag(tB22)^(-1/2),m0,m0)
	corr=diagB22%*%tB22%*%diagB22
	round(corr,3)
}

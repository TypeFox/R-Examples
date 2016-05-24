tb2.pois <-
function(alpha,theta)
{
#This function computes \tilde B_{22} under Poisson mixture.
#\tilde B_{22} has been standardized to the correlation matrix.
#
#alpha:  vector of mixture probabilities.
#theta:  vector of parameters of each component.
	m0=length(alpha)
	size=min(qpois(1-(1e-100),max(theta)),200)
	x=0:size
	pdf=0
	delta=c()
	y=c()
	z=c()

	for(i in 1:m0)
	{
		delta=cbind(delta,dpois(x,theta[i])-dpois(x,theta[m0]))
		y=cbind(y,(x/theta[i]-1)*dpois(x,theta[i]))
		z=cbind(z,(x*(x-1)/theta[i]^2-2*x/theta[i]+1)*dpois(x,theta[i])/2)
		pdf=pdf+alpha[i]*dpois(x,theta[i])
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

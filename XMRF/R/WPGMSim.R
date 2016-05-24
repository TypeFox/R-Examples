WPGMSim <-
function(n, p, R, alpha, Theta, maxit=10000)
{
	X = matrix(rpois(n*p,1),n,p);
	iter = 1;
	while(iter<maxit)
	{
		for(j in 1:p)
		{
			num = exp( matrix(1,n,1)%*%t(alpha[j]*c(0:R)-log(factorial(c(0:R)))) + matrix(c(0:R)%x%X[,-j]%*%Theta[-j,j],n,R+1) );
			Pmat = num/matrix(apply(num,1,sum),n,R+1);
			X[,j] = apply(apply(Pmat,1,mymult)==1,2,which) - 1;
		}
		iter = iter + 1;
	}
	return(X)
}

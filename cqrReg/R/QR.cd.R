QR.cd = function(X, y, tau, beta, maxit, toler)
{

if(missing(toler)){
	toler = 1e-3
	
}

if(missing(maxit))
{
maxit = 200
}


if(missing(beta)){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

betah=QRCDCPP(X,y,beta,toler,maxit,tau)
beta=betah
b=quantile(y-X%*%beta,tau)
return(list(beta=beta,b=b))
}

cqr.mm = function(X, y, tau, beta, maxit, toler)
{

if(missing(maxit)){
	maxit = 200
}

if(missing(toler)){
	toler = 1e-3
	
}
n=dim(X)[1]
k=length(tau)

if(missing(beta)){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

beta=CQRMMCPP(X,y,beta,toler,maxit,tau)
b=quantile(y-X%*%beta,tau)

return(list(beta=beta,b=b))
}

cqr.lasso.mm =  function(X,y,tau,lambda,beta,maxit,toler)
{

if(missing(toler)){
	toler = 1e-3
	
}

if(missing(maxit)){
	maxit = 200
}


n=dim(X)[1]
k=length(tau)

if(missing(beta)){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

if(missing(lambda)){
	lambda=1
}

beta0=CQRMMCPP(X,y,beta,toler,maxit,tau)
beta=CQRPMMCPP(X,y,beta,beta0,toler,maxit,tau,lambda)
b=quantile(y-X%*%beta,tau)
return(list(beta=beta,b=b))
}

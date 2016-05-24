QR.lasso.mm = function(X,y,tau,lambda,beta,maxit,toler)
{

if(missing(toler)){
	toler = 1e-3
	
}

if(missing(maxit)){
maxit = 200
}


n=dim(X)[1]
x=cbind(rep(1,n),X)

if(missing(beta)){
	beta=solve(t(x)%*%x,t(x)%*%y)
}

if(missing(lambda)){
	lambda=1
}

beta0=QRMMCPP(X,y,beta,toler,maxit,tau)[-1]
betah=QRPMMCPP(X,y,beta,beta0,toler,maxit,tau,lambda)
beta=matrix(betah[-1])
b=betah[1]
return(list(beta=beta,b=b))
}

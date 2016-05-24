QR.admm = function(X, y, tau,rho,beta, maxit,toler)
{

if(missing(maxit)){
	maxit = 200
}

if(missing(toler)){
	toler = 1e-3
	
}

n=dim(X)[1]
x=cbind(rep(1,n),X)

if(missing(beta)){
	beta=solve(t(x)%*%x,t(x)%*%y)
}

if(missing(rho)){
	rho=0.4
}

betah=QRADMMCPP(X,y,beta,toler,maxit,tau,rho)
beta=matrix(betah[-1])
b=betah[1]
return(list(beta=beta,b=b))
}

cqr.lasso.admm=  function(X,y,tau,lambda,rho,beta,maxit){
	# # X: n by p design matrix
	# # y: n by 1 response vector
	# # tau: Quantile level
	# # lambda: penalty parameter
	# # rho: augmented lagrangian  parameter
	# # alpha: over-relaxation  parameter
if(missing(maxit)){
	maxit = 200
}

if(missing(rho)){
	rho=0.4
}

if(missing(beta)){
	beta=solve(t(X)%*%X,t(X)%*%y)
}

if(missing(lambda)){
	lambda=1
}

n=dim(X)[1]
p = dim(X)[2]
k=length(tau)

b=mean(y-X%*%beta)
beta=c(rep(b,k),beta)    
	

	bmatrix=matrix(c(rep(1,n),rep (0,n*k)),n*k,k+1)[,-(k+1)]
	X.star= cbind(bmatrix,t(matrix(t (X),p,n*k)))
	y.star=rep(y,k)
	tau.star=c(t(matrix(tau,k,n)))     
      
	betai=CQRPADMMCPP(X.star,y.star,beta,maxit,tau.star,rho,lambda,p,k)
	beta=matrix(betai[-(1:k)])
        b=betai[1:k]
	return(list(beta=beta,b=b))	
}

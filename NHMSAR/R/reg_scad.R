soft_tresh = function(z,lambda) {
	S = 0*z
	S[z>lambda] = z-lambda
	S[z<-lambda] = z+lambda
    return(S)
}
fscad = function(z,lambda,gamma){
	beta = z
	beta[abs(z) <= gamma*lambda] = soft_tresh(z,gamma*lambda/(gamma-1))/1-1/(gamma-1)
	beta[abs(z) <= 2*lambda] = soft_tresh(z,lambda)
	return(beta)
}


reg_scad = function(y,X,lambda,gamma=3.7,weights = 1,beta=0,max.iter=100,eps=1e-4) {
	
	d = dim(X)[2]
	n = dim(X)[1]

	if (length(weights)==1) {weights = rep(weights,n)}
	sw = sum(weights)
	X.sc = X
	s_x = NULL
	for (j in 1:d) {
		mu_x = sum(weights*X[,j])/sw
		s_x[j] = sqrt(sum(weights*X[,j]^2)/sw-mu_x^2)
		X.sc[,j] = (X[,j]-mu_x)/s_x[j]
	}
	X.sc.w = X.sc*matrix(weights,n,d)
	wXX = NULL
	for (j in 1:d) {wXX[j] = t(X.sc.w[,j])%*%X.sc[,j]}
	
	mu = sum(weights*y)/sw
	r = (y-mu)

	if (length(beta)==1){beta = rep(beta,d)	}
	err = 1
	kiter=0
	while (kiter < max.iter & err>eps) { # Si lambda=0 on retrouve la solution de lm
		beta.old = beta
		kiter = kiter+1
		for (j in 1:d) {
			zj = t(X.sc.w[,j])%*%r/wXX[j]+beta[j]
			beta[j] = fscad(zj,lambda,gamma)
			r = r - (beta[j]-beta.old[j])*X.sc[,j]
		}
		err = sum(abs(beta-beta.old))
	}
	
	beta = beta/s_x
	beta0 = mu-sum(beta*apply(X,2,mean))
	return(list(beta0=beta0,beta=beta,y=y,X=X))
}

### ADMM for FGL:
admm.iters = function(Y,lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5,warm=NULL)
{
	K = length(Y)
	p = dim(Y[[1]])[2]
	n=weights

	ns = c(); for(k in 1:K){ns[k] = dim(Y[[k]])[1]}
	S = list(); for(k in 1:K){S[[k]] = cov(Y[[k]])*(ns[k]-1)/ns[k]}

	# initialize theta:
	theta = list()
	for(k in 1:K){theta[[k]] = diag(1/diag(S[[k]]))}
	# initialize Z:
	Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
	# initialize W:
	W = list();	for(k in 1:K) {W[[k]] = matrix(0,p,p) }

	# initialize lambdas:  (shouldn't need to do this if the function is called from the main wrapper function, JGL)
	lam1 = penalty.as.matrix(lambda1,p,penalize.diagonal=penalize.diagonal)
	if(penalty=="fused") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=TRUE)}
	if(penalty=="group") {lam2 = penalty.as.matrix(lambda2,p,penalize.diagonal=penalize.diagonal)}
	# iterations:
	iter=0
	diff_value = 10
	while((iter==0) || (iter<maxiter && diff_value > tol))
	{
		# reporting
#	if(iter%%10==0)
		if(FALSE)
		{
			print(paste("iter=",iter))
			if(penalty=="fused")
			{
				print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
				print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
			}
			if(penalty=="group"){print(paste("crit=",gcrit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))}
		}

		# update theta:
		theta.prev = theta
		for(k in 1:K){
			edecomp = eigen(S[[k]] - rho*Z[[k]]/n[k] + rho*W[[k]]/n[k])
			D = edecomp$values
			V = edecomp$vectors
			D2 = n[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[k]) )
			theta[[k]] = V %*% diag(D2) %*% t(V)
		}

		# update Z:
		# define A matrices:
		A = list()
		for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
		if(penalty=="fused")
		{
			# use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
			if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
			if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}  # the option to not penalize the diagonal is exercised when we initialize the lambda matrices
		}
		if(penalty=="group")
		{
			#  minimize rho/2 ||Z-A||_F^2 + P(Z):
	   		Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
		}

		# update the dual variable W:
		for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

		# bookkeeping:
		iter = iter+1
	  	diff_value = 0
	       for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
		# increment rho by a constant factor:
		rho = rho*rho.increment
	}
	diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
	out = list(theta=theta,Z=Z,diff=diff,iters=iter)
	return(out)
}



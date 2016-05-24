### ADMM for FGL when Theta is known to be diagonal (unconnected nodes):

#	call:	admm.iters.unconnected(Yu,lambda1=lam1.unconnected,lambda2=lam2.unconnected,rho,weights,
#   penalize.diagonal=FALSE,maxiter,tol)
# lam1.unconnected, lam2.u...  vectors of penalties.


admm.iters.unconnected = function(Y,lambda1,lambda2,penalty="fused",rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5,warm=NULL)
{
	K = length(Y)
	for(k in 1:K){Y[[k]]=as.matrix(Y[[k]])}
	p = dim(Y[[1]])[2]
	n=weights
  
 	# define S vectors: empirical covs of each element
	S = list()
	for(k in 1:K)
	{
  	ntemp=dim(Y[[k]])[1]
  	S[[k]]=rep(0,p)
  	for(j in 1:p)
  	{
  		S[[k]][j] = var(Y[[k]][,j])*(ntemp-1)/ntemp
  	}
	}

	# initialize theta:
	theta = list()
	for(k in 1:K){theta[[k]] = 1/S[[k]]}
	# initialize Z:
	Z = list(); for(k in 1:K){Z[[k]] = rep(0,p)}
	# initialize W:
	W = list();	for(k in 1:K){W[[k]] = rep(0,p)}

	# initialize lambdas:
	lam1 = lambda1
	lam2 = lambda2

	# iterations:
	iter=0
	diff_value = 10
	while((iter==0) || (iter<maxiter && diff_value > tol))
	{
	# update theta to minimize -logdet(theta) + <S,theta> + rho/2||theta - Z + W ||^2_2:
	theta.prev = theta
	for(k in 1:K)
	{
		B = n[k]*S[[k]] - rho*(Z[[k]] - W[[k]])  
		theta[[k]] = 1/(2*rho) * ( -B + sqrt(B^2 + 4*rho*n[k]) )  
#		B = S[[k]] - rho/n[k]*(Z[[k]] - W[[k]])  
#		theta[[k]] = n[k]/(2*rho) * ( -B + sqrt(B^2 + 4*rho/n[k]) )  # written as in our paper
	}

	# update Z:
	# define A matrices:
	A = list()
	for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
	if(penalty=="fused")
	{
		# use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
		if(K==2){Z = flsa2(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
		if(K>2){Z = flsa.general(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
	}
	if(penalty=="group")
	{
		 # minimize rho/2 ||Z-A||_F^2 + P(Z):
		Z = dsgl(A,rho,lam1,lam2,penalize.diagonal=TRUE)
	}

	# update the dual variable W:
	for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}
	
	# bookkeeping:
	iter = iter+1
  	diff_value = 0
	for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
	# increment rho by a constant factor:
	rho = rho * rho.increment
	}
	diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
	out = list(theta=theta,Z=Z,diff=diff,iters=iter)
	return(out)
}


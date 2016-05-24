fem.loglik <-
function(prms,K,Y,U){
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	prop = c(prms[1:(K-1)],1-sum(prms[1:(K-1)]))
	alpha = prms[K:(2*K-1)]
	beta = prms[(2*K):(3*K-1)]
	my = matrix(prms[(3*K):length(prms)],ncol=p,byrow=T)
	d = min(p-1,(K-1))

	QQ = matrix(NA,n,K)
	QQ2 = matrix(NA,n,K)
	T = matrix(NA,n,K)
	D = array(0,c(K,p,p))

	# Compute posterior probabilities
	for (k in 1:K){
		D[k,,] = diag(c(rep(alpha[k],d),rep(beta[k],(p-d))))
		bk = D[k,p,p]
		YY = Y-t(matrix(rep(my,n),p,n)) 
		projYY = YY %*% U %*% t(U) # proj dans l'espace latent
		
		if (d==1){
			  for (i in 1:n){QQ[i,k] =  1/D[k,1,1] * sum(projYY[i,]^2) + 1/D[k,p,p]*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)}
		}
		else{
		      sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
		      for (i in 1:n){QQ[i,k] =  projYY[i,] %*% sY %*% as.matrix(projYY[i, ],p,1) + 1/bk*sum((YY[i,] - projYY[i,])^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)}
		}
	}
	A = -1/2 * QQ
 	loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
}


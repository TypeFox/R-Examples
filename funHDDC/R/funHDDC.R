funHDDC <-
function(fd,K,init='hclust',model='AkBkQkDk',thd=0.05,maxit=50,eps=1e-6,...){
# Usage: res <- hddc(Y,K,...){
	# Initialization
	X = t(fd$coefs)
	n = nrow(X)
	p = ncol(X)

	# New objects
	L = rep(c(-Inf),1,(maxit+1))

	# Initialization of T	
  if (init=='kmeans'){
		T = matrix(0,n,K)
		ind = kmeans(X,K)$cluster
		for (i in 1:n){ T[i,ind[i]] = 1 }
	}
	else if (init=='random'){
		T = t(rmultinom(n,1,c(rep(1/K,K))))
	}
	else if (init=='hclust'){
	  T   = matrix(0,n,K)
	  ind = cutree(hclust(dist(X),method='ward'),K)
	  for (k in 1:K){ T[which(ind==k),k] = 1}
	}
	
	prms = list()
	# Main loop
	for (i in 1:maxit){
		cat('*')
		prms.old = prms
		prms = .mstep(fd,T,thd,model)
		T.old = T
		T = .estep(prms,fd)
		L[i+1] = .loglikelihood(prms,fd)
		if (abs(L[i+1]-L[i])<n * eps) break
	}
	cat('\n')

	# Returning the results
	if (L[i+1] < L[i]){cat('Re-use last valid parameters!\n'); prms = prms.old; T = T.old}
	cls = max.col(T)
	crit = .bic(L[(i+1)],prms,n,T);
	res = list(cls=cls,P=T,prms=prms,bic=crit$bic,aic=crit$aic,icl=crit$icl,loglik=L[2:(i+1)])
}

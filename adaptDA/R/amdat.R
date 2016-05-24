amdat <-
function(X,cls,Y,K,model='qda',maxit=25,eps=1e-4,...){
	# Initialization
	C = max(cls)
	p = ncol(X)
	nx = nrow(X)
	ny = nrow(Y)
	n = nx + ny

	# New objects
	L = rep(c(-Inf),1,(maxit+1))
	
	# Case K < C 
	if (K < C) stop('amda: The number of components must be higher or equal to ',C,'!')

	# Case K >= C 
	if (K >= C){
		#cat('amda: Classification with unobserved classes (K>=C)\n')
		# Initialization
		prms = .amdat.initiparam(X,cls,Y,K)
		S = matrix(0,nx,K)
		for (i in 1:K) {S[cls==i,i] = 1}
		# Main loop
		for (i in 1:maxit){
			#cat('.')
			P = .amdat.estep2(prms,X,S,Y,C,K)
			prms = .amdat.mstep(P,X,Y,C,K)
			if (sum(prms$prop < 1e-3) == 1){
				cat('!!! Empty class !!!\n')
				L[i+1] = -Inf
				break
			}
			L[i+1] = .amdat.loglikelihood(prms,P,X,Y,C,K)
			if (abs(L[i+1]-L[i])<n*eps) break
		}
		#cat('\n')
	}
	
	# Returning the results
  Py = P[(nx+1):n,]
	cls = max.col(Py) # Prediction for new data
	crit = .amdat.crit(L[(i+1)],P,prms,n);
	res = list(model=model,K=prms$K,p=p,mean=prms$mean,prop=prms$prop,V=prms$V,P=P,Py=Py,cls=cls,crit=crit,loglik=L[2:(i+1)])
	class(res) <- "amdat"
	res
}

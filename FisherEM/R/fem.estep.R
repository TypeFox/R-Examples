fem.estep <-
function(prms,Y,U){
# require('MASS')
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	K = prms$K
	mu = prms$mean
	prop = prms$prop
	D = prms$D
	d = min(p-1,(K-1))
	QQ = matrix(NA,n,K)
	T = matrix(NA,n,K)
	memlim = 2500
	
	# Compute posterior probabilities
	for (k in 1:K){
		bk = D[k,p,p]
		mY = prms$my[k,]
		YY = Y-t(matrix(rep(mY,n),p,n)) 
		projYY = YY %*% U %*% t(U)
		if (nrow(YY)<=memlim){
			if (d==1){
				QQ[,k] =  1/D[k,1,1] * rowSums(projYY^2) + 1/D[k,p,p]*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)
			} else{
				sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
				QQ[,k] =   diag(projYY %*% sY %*% t(projYY)) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)
			}
		}
		else{ # Do the same but by blocks of memlim obs.
			if (n%%memlim==0) quot = n %/% memlim # quotient de la division euclidienne
      else quot = n %/% memlim + 1
			for (i in 1:quot){
				if (d==1){
					QQ[(memlim*(i-1)+1):min(memlim*i,n),k] =  1/D[k,1,1] * rowSums(projYY[(memlim*(i-1)+1):min(memlim*i,n),]^2) + 1/D[k,p,p]*rowSums((YY[(memlim*(i-1)+1):min(memlim*i,n),] - projYY[(memlim*(i-1)+1):min(memlim*i,n),])^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)
				} else{
					sY = U %*% ginv(D[k,(1:d),(1:d)]) %*% t(U)
					QQ[(memlim*(i-1)+1):min(memlim*i,n),k] =   diag(projYY[(memlim*(i-1)+1):min(memlim*i,n),] %*% sY %*% t(projYY[(memlim*(i-1)+1):min(memlim*i,n),])) + 1/bk*rowSums((YY[(memlim*(i-1)+1):min(memlim*i,n),] - projYY[(memlim*(i-1)+1):min(memlim*i,n),])^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)
				}
			}
		
		}
	}
	# Compute the log-likelihood
	A = -1/2 * QQ
 	loglik = sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))
	for (k in 1:K) {T[,k] = 1 / rowSums(exp(0.5*(QQ[,k]*matrix(1,n,K)-QQ)))}
	
	# Return the results
	list(T=T,loglik=loglik)
}


fem.mstep <-
function(Y,U,T,model,method){
	# 12 different submodels: [DkBk] ... [AkjBk]
	# Initialization
	Y = as.matrix(Y)
	n = nrow(Y)
	p = ncol(Y)
	K = ncol(T)
	d = min(p-1,(K-1))
	U = as.matrix(U)

	mu   = matrix(NA,K,d)
	m   = matrix(NA,K,p)
	prop = rep(c(NA),1,K)
	D = array(0,c(K,p,p))

	# Projection
	X = as.matrix(Y) %*% as.matrix(U)

	# Estimation
	test = 0
	for (k in 1:K){
		
		nk  = sum(T[,k])
		if (nk ==0) stop("some classes become empty\n",call.=FALSE)
		# Prior Probability
		prop[k] = nk / n
		# Mean in the latent space
		mu[k,]  = colSums((T[,k]*matrix(1,n,d))* X) / nk
		# Observed space
		m[k,]  = colSums(T[,k]*matrix(1,n,p)* Y) / nk
		YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
		if (nk < 1) denk = 1 else denk = (nk-1)
		Ck  = crossprod(T[,k]*matrix(1,n,p)* YY, YY) / denk #crossprod(x,y) = t(x) %*% y
		C   = cov(Y)

		# Estimation of Delta k amongst 8 submodels
		if (model=='DkBk'){
			D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		if (model=='DkB'){
			D[k,(1:d),(1:d)] = crossprod(Ck%*%U,U)
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		if (model=='DBk'){
			D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		if (model=='DB'){
			D[k,(1:d),(1:d)] = crossprod(C%*%U,U)
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		
		if (model=='AkjBk'){
		if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		
		if (model=='AkjB'){
		if (d==1){D[k,1,1] = diag(crossprod(Ck%*%U,U))} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(Ck%*%U,U)))}
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		
		if (model=='AkBk'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		
		if (model=='AkB'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(Ck%*%U,U)))/d} else{
			  D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(Ck%*%U,U)))/d,d))}
			  bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			  bk[bk<=0] = 1e-3
			  if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		
		if (model=='AjBk'){
		if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else {
			D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		if (model=='AjB'){
		if (d==1){D[k,1,1] = diag(crossprod(C%*%U,U))} else{
			D[k,(1:d),(1:d)] = diag(diag(crossprod(C%*%U,U)))}
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		if (model=='ABk'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
			bk = (sum(diag(Ck)) - sum(diag(crossprod(Ck%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
		if (model=='AB'){
		if (d==1){D[k,1,1] = sum(diag(crossprod(C%*%U,U)))} else {
			D[k,(1:d),(1:d)] = diag(rep(sum(diag(crossprod(C%*%U,U)))/d,d))}
			bk = (sum(diag(C)) - sum(diag(crossprod(C%*%U,U)))) / (p-d)
			bk[bk<=0] = 1e-3
			if ((d+1)==p) D[k,p,p] = bk else D[k,((d+1):p),((d+1):p)] = diag(rep(bk,p-d))	
			if (sum(diag(D[k,,]<1e-3))!=0) test = test+1
		}
	}
	prms = list(K=K,p=p,mean=mu,my=m,prop=prop,D=D,model=model,method=method,V=U,test=test)
}


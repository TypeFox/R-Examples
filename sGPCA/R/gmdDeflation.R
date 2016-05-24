gmdDeflation <- function(X,Q,R,k,n,p){
	
##computation
	
	ugmd = matrix(nrow = n,ncol = k)
	vgmd = matrix(nrow = p, ncol = k)
	dgmd = rep(0,k)
	propv = rep(0,k)
	Xhat = X
	
	u = rnorm(n)
	v = rnorm(p)
	
	qrnorm = sum(diag(t(X) %*% Q %*% X %*% R))
	
	cumv = rep(0,k)
	
	thr = 1e-6
	for(i in 1:k){
		err=1
		while(err > thr){
			oldu = u
			oldv = v
			uhat = Xhat %*% R %*% v
			u = uhat/as.double(sqrt(t(uhat)%*% Q %*% uhat))
			vhat = t(Xhat) %*% Q %*% u
			v = vhat/as.double(sqrt(t(vhat) %*% R %*% vhat))
			err = t(oldu - u) %*% (oldu - u) + t(oldv -v ) %*% (oldv - v)
		}
		dgmd[i] = t(u) %*% Q %*% X %*% R %*% v
		ugmd[,i] = u
		vgmd[,i] = v
		Xhat = Xhat - dgmd[i] *  u %*% t(v)
		propv[i] = dgmd[i]^2/as.double(qrnorm)
		cumv[i] = sum(propv[1:i])
	}
	
	
	return.object = list(ugmd,vgmd[,1:k],dgmd,cumv,propv)
	return(return.object)
	
}



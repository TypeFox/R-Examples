
hdmda <- function(X,cls,K=1:10,model='AkjBkQkDk',show=FALSE,...){
	cls = factor(cls)
	if (!is.numeric(cls)) {
		levels = levels(cls)
		cls = as.numeric(cls)
	} else levels = NULL
	C = max(cls)
	alpha = rep(NA,C); prms = list()
	for (c in 1:C){
		alpha[c] = sum(cls==c) / nrow(X)
		prms[[c]] = hddc(X[cls==c,],K=K,model=model,show=FALSE,...)
	}
	obj = list(alpha=alpha,prms=prms,kname=levels)
	class(obj) = 'hdmda'
	return(obj)
}

predict.hdmda <- function(object,X,...){
	p = ncol(X)
	N = nrow(X)
	C = length(object$alpha)
	P = matrix(NA,nrow(X),C)
	for (c in 1:C){
		par = object$prms[[c]]
		K <- par$K
		a <- par$a
		b <- par$b
		mu <- par$mu
		d <- par$d
		prop <- par$prop
		Q <- par$Q
		
		b[b<1e-6] <- 1e-6
		
		if(par$model=="AJBQD") {
			K_pen <- diag((mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(mu)))-2*(mu%*%Q%*%diag(1/a[1,1:d[1]],d[1]))%*%(t(Q)%*%t(X))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(X)+2*(mu%*%Q)%*%(t(Q)%*%t(X))-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))
		}
		else if(par$model=="ABQD") {
			K_pen <- diag(1/a[1]*(mu%*%Q)%*%(t(Q)%*%t(mu)))+1/b[1]*(diag(tcrossprod(mu))-2*mu%*%t(X)-diag(tcrossprod(mu%*%Q)))-2*log(c(prop))+2*(1/b[1]-1/a[1])*(mu%*%Q)%*%(t(Q)%*%t(X))
		}
		else{
			K_pen <- matrix(0,K,N)
			for (i in 1:K) {
				s <- sum(log(a[i,1:d[i]]))
				Xk <- as.matrix(X-matrix(mu[i,],N,p,byrow=TRUE))
				proj <- (Xk%*%Q[[i]])%*%t(Q[[i]])
				A <- (-proj)%*%Q[[i]]%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
				B <- Xk-proj
				K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
			}
		}
		K_pen = -1/2*t(K_pen)
		P[,c] = object$alpha[c] * rowSums(exp(K_pen))
	}
	PP <- matrix(NA,N,C)
	for (c in 1:C) {PP[,c] = P[,c] / rowSums(P)}
	if (!is.null(object$kname)) {class = object$kname[max.col(PP)]}
	else {class = max.col(PP)}
	list(class=class,posterior=PP)
}


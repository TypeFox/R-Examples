##
off <- function(V,C) {
	F <- V %*% C %*% t(V)
	sum(diag(t(F) %*% F)) - sum(diag(F^2))
}

ffdiag <- function(C0, V0=NULL,eps=.Machine$double.eps,
		itermax=200, keepTrace=FALSE) {
	n <- dim(C0)[1]
	K <- dim(C0)[3]
	C = C0
	Id=diag(n)
	V <- V0
	if (is.null(V)) V <- Id
	inum = 0
	df = 1
	crit <- NA
	B_trace <- V
	while ((df > eps) & (inum < itermax)) {
		W = matrix(.C("getW",
			as.double(C),
			as.integer(n),
			as.integer(K),
			as.double(diag(0,n,n)))[[4]],n,n,byrow=TRUE)
		e <- ceiling(log2(max(colSums(abs(t(W))))))
		s = max(0,e-1)
		W = W/(2^s)
		V = (Id+W) %*% V
		V = diag(1/sqrt(diag(V %*% t(V)))) %*% V
		C <- array(unlist(sapply(1:K, function(i) V %*% C0[,,i] %*%
			t(V),simplify=F)),dim=c(n,n,K))
		f <- sum(apply(C0, 3, function(X) off(V,X)))
		crit <- c(crit,f)
		if (inum > 2) df = abs(crit[inum-1]-crit[inum])
		inum <- inum +1
		if(keepTrace) B_trace <- cbind(B_trace,V)
	}
	if (inum == itermax) {
		warning("Convergence not reached")
	}
	if (keepTrace) {
		return(list(B=V,criter=crit,
			B_trace=array(B_trace,dim=c(n,n,inum+1))))
	}
	else {
		return(list(B=V,criter=crit,B_trace=NULL))
	}
}

			
	
	

































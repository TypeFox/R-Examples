## 
uwedge <- function(M,W_est0=NULL,eps=.Machine$double.eps,
		itermax=200,keepTrace=FALSE) {
	d <- dim(M)[1]
	L <- dim(M)[3]
	iter <- 0
	improve <- 10
	W_est <- W_est0
	if (is.null(W_est)) {
		eigendec <- eigen(M[,,1])
		W_est <- diag(1/sqrt(abs(eigendec$values))) %*% t(eigendec$vectors)
	}
	B_trace <- W_est
	M <- array(unlist(sapply(1:L, function(l)
					0.5*(M[,,l]+t(M[,,l])),simplify=F)),dim=c(d,d,L))		
	Ms <- array(unlist(sapply(1:L, function(l)
					W_est %*% M[,,l] %*% t(W_est),simplify=F)),dim=c(d,d,L))
	Rs <- sapply(1:L, function(l) diag(Ms[,,l]))
	crit <- sum(Ms^2)-sum(Rs^2)
	while (improve>eps & iter<itermax) {
		B=Rs %*% t(Rs)
		C1 <- sapply(1:d, function(daux) rowSums(Ms[,daux,] * Rs))
		D0= B*t(B)-diag(B) %*% t(diag(B))
		A0=diag(1,d)+(C1 * B - diag(diag(B))%*% t(C1))/(D0+diag(1,d))
		W_est=solve(A0) %*% W_est
		Raux=W_est %*% M[,,1] %*% t(W_est)
		aux=1/sqrt(abs(diag(Raux)))
		W_est=diag(aux) %*% W_est
		Ms <- array(unlist(sapply(1:L, function(l)
					W_est %*% M[,,l] %*% t(W_est),simplify=F)),dim=c(d,d,L))
		Rs <- sapply(1:L, function(l) diag(Ms[,,l]))
		critic <- sum(Ms^2)-sum(Rs^2)
		improve <- abs(critic-rev(crit)[1])
		crit <- c(crit,critic)
		iter <- iter+1
		if(keepTrace) B_trace <- cbind(B_trace,W_est)
	}
	if (iter == itermax) {
		warning("Convergence not reached")
	}
	if (keepTrace) {
		return(list(B=W_est,criter=crit,
			B_trace=array(B_trace,dim=c(d,d,iter+1))))
	}
	else {
		return(list(B=W_est,criter=crit,B_trace=NULL))
	}
}














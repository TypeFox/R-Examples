err <- function(W,C0,C,p) {
	N <- dim(C)[1]
	M <- N
	T <- dim(C)[3]
	E = 0;
	W = W %*% diag(1/sqrt(diag(t(W) %*% C0 %*% W)))
	for (t in 1:T) {
		D = t(W) %*% C[,,t] %*% W
	    D = D-diag(diag(D));
		E = E+p[t]*sum(D^2)
	}
	E = E/(M^2-M)
}

qdiag <- function(C, W0=NULL, eps=.Machine$double.eps,
		itermax=200, keepTrace=FALSE) {
	N <- dim(C)[1]
	T <- dim(C)[3]
	W <- W0
	if (is.null(W)) W <- matrix(rnorm(N^2),N,N)
	p = rep(1/T,T)
	C0 <- diag(N)
	crit <- NULL
	deltacrit <- 1
	W = W %*% diag(1/sqrt(diag(t(W) %*% C0 %*% W)))
	eigC0 <- eigen(C0)
	P <- diag(1/sqrt(eigC0$values)) %*% t(eigC0$vectors)
	C <- array(unlist(sapply(1:T, function(i) P %*% C[,,i] %*%
		t(P),simplify=F)),dim=c(N,N,T))
	issymmetric <- apply(C,3,function(X) isSymmetric(X))
	C0 = P %*% C0 %*% t(P)
	W  = t(solve(P)) %*% W
	# initialisations for OKN^3
	D <- matrix(0,N,N)
	for (t in 1:T) {
		M1 = C[,,t] %*% W
		if (issymmetric[t])	D  = D + (2*p[t]*(M1 %*% t(M1)))
		else {
			M2 = t(C[,,t]) %*% W
			D  = D + (p[t]*( M1%*%t(M1) + M2%*%t(M2) ))
		}
	}
	n <- 0
	B_trace <- W
	while (n < itermax & deltacrit > eps) {
		delta_w = 0
		for (i in 1:N) {
			w = matrix(W[,i],N,1)
			for (t in 1:T) {
				m1 = C[,,t] %*% w
				if (issymmetric[t])	D  = D - (2*p[t]*(m1 %*% t(m1)))
				else {
					m2 = t(C[,,t]) %*% w
					D  = D - (p[t]*( m1%*%t(m1) + m2%*%t(m2) ))
				}
			}
			eigD <- eigen(D)
			w_new = eigD$vectors[,N]
			delta_w = max( delta_w, min(sqrt(sum((w-w_new)^2)),
						sqrt(sum((w+w_new)^2))))
			for (t in 1:T) {
				m1 = C[,,t] %*% w_new
				if (issymmetric[t])	D  = D + (2*p[t]*(m1 %*% t(m1)))
				else {
					m2 = t(C[,,t]) %*% w_new
					D  = D + (p[t]*( m1%*%t(m1) + m2%*%t(m2) ))
				}
			}
			W[,i] <- w_new
		}
		crit <- c(crit,err(W,C0,C,p))
		if (n>1) deltacrit <- abs(crit[n]-crit[n-1])
		n <- n+1
		if (keepTrace) B_trace <- cbind(B_trace,t(t(P) %*% W))
	}
	W = t(t(P) %*% W)
	if (n == itermax) {
		warning("Convergence not reached")
	}
	if (keepTrace) {
		return(list(B=W,criter=crit,
			B_trace=array(B_trace,dim=c(N,N,n+1))))
	}
	else {
		return(list(B=W,criter=crit,B_trace=NULL))
	}
}

		





























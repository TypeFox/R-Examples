# Auxiliary function:

SqrtS <- function(S)
# For a symmetric, positive definite matrix S,
# this procedure returns a matrix B such that
#    S = B %*% t(B).
{
	res <- eigen(S,symmetric=TRUE)
	B <- t(t(res$vectors) * sqrt(res$values))
	return(B)
}

# Function to check how well a matrix parameter
# satisfies the fixed-point equation:

CheckS0 <- function(X,nu,S)
{
	n <- dim(X)[1]
	p <- dim(X)[2]
	denom <- nu + rowSums(X*t(qr.solve(S,t(X))))
	Z <- X/sqrt(denom)
	S1 <- (nu + p) * crossprod(Z)/n
	err1 <- norm(S-S1)
	R <- SqrtS(S)
	Y <- t(qr.solve(R,t(X)))
	denom <- nu + rowSums(Y^2)
	Z <- Y/sqrt(denom)
	S2 <- (nu + p) * crossprod(Z)/n
	err2 <- norm(S2 - diag(rep(1,p)))
	return(c(err1=err1,err2=err2))
}


# Algorithm for the location-scatter problem:

MVTMLEr <- function(X,nu=1,delta=10^(-7))
{
	if (nu < 1)
	{
		print("Need nu >= 1 degress of freedom!")
		print("Use nu == 1 now!")
		nu <- 1
	}
	Xnames <- colnames(X)
	X <- as.matrix(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	Y <- cbind(X,rep(1,n))
	G <- MVTMLE0r(Y,nu-1,delta)$S
	G <- G / G[p+1,p+1]
	mu.hat <- G[1:p,p+1]
	names(mu.hat) <- Xnames
	Sigma.hat <- G[1:p,1:p] - tcrossprod(mu.hat)
	rownames(Sigma.hat) <-
		colnames(Sigma.hat) <- Xnames
	Shape.hat <- Sigma.hat / det(Sigma.hat)^(1/p)
	return(list(mu.hat=mu.hat,
		Sigma.hat=Sigma.hat,
		Shape.hat=Shape.hat))
}





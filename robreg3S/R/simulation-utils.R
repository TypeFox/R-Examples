generate.randbeta <- function(p){
	tmp <- matrix(rnorm(p), ncol = p)
	.norm.vec <- function (x) x/sqrt(sum(x^2))
	t(apply(tmp, 1, .norm.vec))
}

generate.cellcontam.regress <- function(n, p, A, sigma, b, k, cp){
	if( cp >= 0.1 ) stop("Too much cellwise outliers!")

	x <- mvrnorm(n, mu=rep(0, p), Sigma=A)
	e <- rnorm(n, sd=sigma)
	y <- c(matrix(b, nrow=1) %*% t(x) + e)

	y[1:round(n*cp)] <- c(matrix(b, nrow=1) %*% t(x[1:round(n*cp),]) + k * sigma )
	cX <- matrix( k, n, p) 
	cU <- rbinom( n*p, 1, cp)
	x <- x * (1 - cU) + cX * cU
	list(y=y, x=x)
}

generate.casecontam.regress <- function(n, p, A, sigma, b, l, k, cp){
	if( cp >= 0.3 ) stop("Too much casewise outliers!")

	x <- mvrnorm(n, mu=rep(0, p), Sigma=A)
	e <- rnorm(n, sd=sigma)
	y <- c(matrix(b, nrow=1) %*% t(x) + e)

	Aeig <- eigen(A, symmetric=T)$vector
	Aevec <- Aeig[,p]
	Aevec.size <- sqrt(Aevec%*%solve(A)%*%Aevec)
	Aevec <- l*Aevec/Aevec.size
	if(cp > 0){ 
		x[1:round(cp*n),] <- matrix( Aevec, round(cp*n), p, byrow=TRUE) 
		y[1:round(cp*n)] <- c(matrix(b, nrow=1) %*% t(x[1:round(cp*n),]) ) + rnorm(round(cp*n), mean=k, sd=sigma)
	}
	list(y=y, x=x)
}


generate.cellcontam.regress.dummies <- function(n, p, pd, probd, A, sigma, b, k, cp){
	if( cp >= 0.1 ) stop("Too much cellwise outliers!")

	x <- mvrnorm(n, mu=rep(0, p + pd), Sigma=A)
	dummies <- x[,(p+1):(p+pd),drop=FALSE]
	x <- x[,1:p]
	for(j in 1:pd) dummies[,j] <- ifelse( dummies[,j] < qnorm( probd[j]), 1, 0)
	e <- rnorm(n, sd=sigma)
	y <- c(matrix(b[1:p], nrow=1) %*% t(x) + matrix(b[(p+1):(p+pd)], nrow=1) %*% t(dummies) + e)

	y[1:round(n*cp)] <- c(matrix(b[1:p], nrow=1) %*% t(x[1:round(n*cp),]) + 
		matrix(b[(p+1):(p+pd)], nrow=1) %*% t(dummies[1:round(n*cp),]) + k * sigma )
	cX <- matrix( k, n, p) 
	cU <- rbinom( n*p, 1, cp)
	x <- x * (1 - cU) + cX * cU
	list(y=y, x=x, dummies=dummies)
}

generate.casecontam.regress.dummies <- function(n, p, pd, probd, A, sigma, b, l, k, cp){
	if( cp >= 0.3 ) stop("Too much casewise outliers!")

	x <- mvrnorm(n, mu=rep(0, p + pd), Sigma=A)
	dummies <- x[,(p+1):(p+pd),drop=FALSE]
	x <- x[,1:p]
	for(j in 1:pd) dummies[,j] <- ifelse( dummies[,j] < qnorm( probd[j]), 1, 0)
	e <- rnorm(n, sd=sigma)
	y <- c(matrix(b[1:p], nrow=1) %*% t(x) + matrix(b[(p+1):(p+pd)], nrow=1) %*% t(dummies) + e)

	A <- A[1:p,1:p]
	Aeig <- eigen(A, symmetric=T)$vector
	Aevec <- Aeig[,p]
	Aevec.size <- sqrt(Aevec%*%solve(A)%*%Aevec)
	Aevec <- l*Aevec/Aevec.size
	if(cp > 0) x[1:round(cp*n),] <- matrix( Aevec, round(cp*n), p, byrow=TRUE) 
	if(cp > 0) y[1:round(cp*n)] <- c(matrix(b[1:p], nrow=1) %*% t(x[1:round(n*cp),]) + 
		matrix(b[(p+1):(p+pd)], nrow=1) %*% t(dummies[1:round(n*cp),]) + rnorm( round(cp*n), k, sigma) )
	list(y=y, x=x, dummies=dummies)
}


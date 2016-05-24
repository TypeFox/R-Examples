
generate.randcorr <- function(cond, p, tol=1e-5, maxits=100) {
	if( p < 3 ) stop("p must be larger than 2.")
	lambda <- sort(c(1, runif(p-2, min=1, max=cond), cond))
	x <- matrix(rnorm(p*p),p,p)
	Sigma <- x%*%t(x)
	Q <- eigen(Sigma, symmetric=TRUE)$vectors
	ratio <- 0
	iter <- 0
	while (abs(ratio - cond) > tol & iter < maxits) {
		iter <- iter + 1
		Sigma <- Q%*%diag(lambda)%*%t(Q)
		Sigma <- diag(diag(Sigma)^(-1/2))%*%Sigma%*%diag(diag(Sigma)^(-1/2))
		eS <- eigen(Sigma, symmetric=TRUE)
		Q <- eS$vectors
		lambda <- eS$values
		ratio <- lambda[1]/lambda[p]
		lambda[p] <- lambda[1]/cond
	}
	return(Sigma)
}

.generate.clean <- function(n, p, cond, A=NULL){
	if( is.null(A) ) A <- generate.randcorr(cond, p)
	x <- mvrnorm(n, mu=rep(0, p), Sigma=A)
	return(list(x=x, A=A))
}

generate.cellcontam <- function(n, p, cond, contam.size, contam.prop, A=NULL){
	x <- .generate.clean(n, p, cond, A)
	contam.num <- floor(n*p*contam.prop)
	u <- matrix( 0, n, p)
	if( contam.num > 0){
		u[sample(1:(n*p), contam.num)] <- 1
		x$x[ which(u == 1)] <- contam.size + rnorm(contam.num, sd=0.01)
	}
	x$u <- u
	x
}

generate.casecontam <- function(n, p, cond, contam.size, contam.prop, A=NULL){
	x <- .generate.clean(n, p, cond, A)
	Aeig <- eigen(x$A, symmetric=T)$vector
	Aevec <- Aeig[,p]
	Aevec.size <- sqrt(Aevec%*%solve(x$A)%*%Aevec)
	Aevec <- contam.size*Aevec/Aevec.size
	contam.num <- floor(n*contam.prop)
	u <- matrix(0, n, p)
	if( contam.num > 0){
		##u[ sample(1:n, contam.num), ] <- 1
		u[ 1:contam.num, ] <- 1
		x$x[ rowSums(u) == p] <- matrix(Aevec, contam.num, p, byrow=T) + rnorm(contam.num*p, sd=0.01)
	}
	x$u <- u
	x
}



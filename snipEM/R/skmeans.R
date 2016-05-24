# compile C code with: 
# R CMD SHLIB sniptotal.c -L. -lvmr 
# then: source("snipC.r")

#dyn.load("sniptotal.dll")

# x: a numerical data matrix (n by d) 
# k: number of clusters, k>=1
# itermax: number of iterations of the algorithm
# V: binary matrix for snipping. starting solution. number of zeros
# will be preserved and correspond to snipped entries. 
# s: binary vector of size n for trimming. starting solution. Number
# of zeros will be preserved and correspond to trimmed rows. NULL performs no trimming
# clust: vector of size n containing values from 1 to k. starting
# solution for class labels. 
# D: tuning parameter for the fitting algorithm. Corresponds
# approximately to the maximal change in loss by switching two non
# outlying entries. Comparing different choices is recommended.  

## WARNING: x and V are matrices of the same size. Not data frames. 

skmeans=function(X, k, V, clust, s, itersmax=10^5, D=1e-1) {
	## check dat
	if( missing(X) ) stop("'X' missing")
	if( missing(V) ) stop("'V' missing")
	if( missing(k) ) stop("'k' missing")
	if( missing(clust) ) stop("'clust' missing")
	
	if(is.data.frame(X) | is.matrix(X))
		X <- data.matrix(X)
	else stop("Data matrix must be of class matrix or data.frame")
	n <- nrow(X)
	d <- ncol(X)
	
	if(missing(s)) {
		s <- rep(1,n)
	}

	estmea <- function(X,s,V,clust) {
		Y <- X
		Y[V==0] <- NA
		apply(Y[s==1,],2,function(x) tapply(x,clust[s==1],mean,na.rm=T))
	}

  mu <- as.matrix(estmea(X,s,V,clust))

	loss <- function(X,mu,s,V,clust) {
		Y <- (X-mu[clust,])^2
		Y[V==0] <- NA
		sum(Y[s==1,],na.rm=T)
	}

	X <- as.vector(t(X))
	V <- as.vector(t(V))
	mu <- as.vector(t(mu))

	out <- .C("sniptotal", x=as.double(X), k=as.integer(k),
		itermax=as.integer(itersmax), V=as.integer(V), s=as.integer(s), 
		clust=as.integer(clust), D=as.double(D), n=as.integer(n),
		d=as.integer(d), mu=as.double(mu))

	V <- matrix(out$V,n,d,byrow=T)

	mu <- estmea( t(matrix(X,d,n)),out$s,V,out$clust)

	return(list(loss=loss(X,mu,out$s,V,out$clust), mu=mu, s=out$s, V=V, clust=out$clust))
}

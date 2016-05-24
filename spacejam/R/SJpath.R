SJpath <- function(X,bfun=function(x){cbind(x,x^2,x^3)}, lambda=NULL,length =10, verbose=FALSE,b0=NULL, maxit=100, tol = .Machine$double.eps^0.25, G.max=NULL){
	cl <- match.call()
	bfun <- match.fun(bfun)
	n <- nrow(X)
	p <- ncol(X)
	k <- ncol(bfun(X[,1]))
	
	if(!is.null(G.max)){
		if(class(G.max) =="igraph") G.max <- get.adjacency(G.max, sparse=TRUE) ==1
		stopifnot(storage.mode(G.max)=="logical")	
	}
	
	#make orthonormal basis functions, and order properly:
	X <- apply(X,2,scale)
	bigX <- matrix(NA,nrow=k*p, ncol=n)
	for(j in 1:p){
		bigX[p*(1:k -1) + j,] <- t(svd(apply(bfun(X[,j]),2,scale,scale=F))$u*sqrt(n-1))
	}
	
	if(is.null(lambdaseq)){
		M <- matrix(NA, nrow=p,ncol=p)
		tmp <- bigX%*%X
		
		for(i in 1:p){
			for(j in (1:p)[-i]){
				M[i,j] <- sqrt(sum(tmp[p*(1:k -1) +j,i]^2))/(n-1)
			}
		}
		M <- M + t(M)
		lmin <- sort(c(M))[2*p]
		lmax <- max(M,na.rm=T)
		lambdaseq <- rev(exp(seq(log(lmin),log(lmax), length = length)))
	}
	
	if(verbose) cat("Progress: \n")
	#fit model
	
	adjmat <- array(dim=c(p,p,length(lambdaseq)))
	for(i in 1:length(lambda)){
		out  <- grpsel(X=X, bigX = bigX, lambda=lambda[i], n=n, k=k, p=p, tol = tol, maxit=maxit,b0=b0, G.max = G.max)
		b0 = out$b0
		adjmat[,,i] <- out$G 
		if(verbose) cat(".")
	}
	#return output
	class(out) <- "SJ"
	out$X = X
	out$graph <- apply(adjmat,3,function(GG){graph.adjacency(GG,mode="undirected")})
	out$cl <- cl
	out$bigX = bigX
	out$n = n
	out$k = k
	out$p = p
	out$tol = tol
	out$maxit = maxit
	out$verbose= verbose
	out$G.max = G.max
	out$lambda =lambda
	out$G = adjmat
	return(out)
}

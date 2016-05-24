SJ <- function(X,bfun=bs, lambda=NULL,length =NULL, verbose=FALSE,b0=NULL, maxit=100, tol = .Machine$double.eps^0.25, G.max=NULL){
	cl <- match.call()
	bfun <- match.fun(bfun)
	n <- nrow(X)
	p <- ncol(X)
	k <- ncol(bfun(X[,1]))
	
	if(!is.null(b0) & any(dim(b0) != c(n,p*k))){
		stop("dimensions of b0 not compatible with basis function")
	}
	
	if(!is.null(G.max)){
		if(class(G.max) =="igraph") G.max <- get.adjacency(G.max, sparse=TRUE) ==1
		stopifnot(storage.mode(G.max)=="logical")	
		stopifnot(ncol(X) == ncol(G.max))
	} 
	
	if(!is.null(lambda) & !is.null(length)){
		if(length(lambda) != length) warning("Lambda and length not compatible. Using lambda")
	}
	
	if(is.null(lambda) & is.null(length)) length = 10
	if(!is.null(lambda) & any(lambda < 0)) stop("lambda must be positive")
	#make orthonormal basis functions, and order properly:
	X <- apply(X,2,scale)
	bigX <- matrix(NA,nrow=k*p, ncol=n)
	
	for(j in 1:p){
		bigX[p*(1:k -1) + j,] <- t(svd(apply(bfun(X[,j]),2,scale,scale=F))$u*sqrt(n-1))
	}
	
	if(is.null(lambda)){
		M <- matrix(NA, nrow=p,ncol=p)
		tmp <- bigX%*%X/(n-1)
		
		for(i in 1:p){
			for(j in (1:p)[-i]){
				M[i,j] <- sum(tmp[p*(1:k -1) +j,i]^2)
			}
		}
		M <- sqrt(M + t(M))
		lmin <- sort(c(M))[p]
		lmax <- max(M, na.rm=T)
		lambda <- rev(log(seq(exp(lmin),exp(lmax), length = length)))
	}
	
	if(verbose) cat("Progress: \n")
	#fit model
	dfs <- rss <- matrix(0,p,length(lambda))	
	adjmat <- array(dim=c(p,p,length(lambda)))
	for(i in 1:length(lambda)){
		out  <- grpsel(X=X, bigX = bigX, lambda=lambda[i], n=n, k=k, p=p, tol = tol, maxit=maxit,b0=b0, G.max = G.max)
		b0 = out$b0
		adjmat[,,i] <- out$G
		tmp <- rowsum(t(out$betas^2)*(n-1), rep(1:p,k))
		dfs[,i] <- colSums(out$G) + (k-1)*rowSums((tmp)/(tmp+lambda[i]))
		rss[,i] <- colSums((X - t(out$betas%*%bigX))^2)
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
	out$dfs = dfs
	out$rss = rss
	out$ecount <- unlist(lapply(out$graph,ecount))
	out$bic = colSums(n*log(rss) + dfs*log(n))
	return(out)
}

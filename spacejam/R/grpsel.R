grpsel <-
function(X, bigX, lambda, n, k, p, tol = .Machine$double.eps^0.25, maxit=100,b0=NULL,verbose=FALSE, G.max = NULL){
	require("Matrix")
	if(is.null(b0)){
		b0 <- matrix(0,nrow=p,ncol=k*p)
		res <- X
	} else {
			res <- X-t(b0%*%bigX)
	}
	if(is.null(G.max)){
		G.old <- matrix(TRUE,nrow=p,ncol=p)
		diag(G.old) <- FALSE
		all.Active <- AtoL(G.old)
		q <- q.max <- dim(all.Active)[2]
	} else {
		G.old <- matrix(TRUE,nrow=p,ncol=p)
		diag(G.old) <- FALSE
		all.Active <- AtoL(G.max)
		q <- q.max <- dim(all.Active)[2]
	}

	if(verbose){
		ssres <- sum(res^2)
		b2 <- b0^2
		pen <- Matrix(0,p,p)
		for(ipen in p*(1:k-1)) pen <- pen + b2[,ipen+1:p] + t(b2[,ipen+1:p])
		pen <- sum(sqrt(pen))/2
		cat("obj0 = ", lambda*pen*(n-1)+(0.5)*ssres)
	}

	its <- 0
	d <- Inf
	fit <- .C("grpsel", b0=as.double(t(b0)), bigX = as.double(bigX),ssx=as.double(1), lambda=as.double(lambda), res=as.double(t(res)),  
    	Active=as.integer(t(all.Active)), n= as.integer(n), p=as.integer(p), k=as.integer(k), q=as.integer(q.max))
	
	b0 <- fit$b0
	G  <- matrix(b0 ,nrow=p,ncol=k*p,byrow=TRUE)[,1:p] != 0
	Active <- AtoL(G)
	q <- length(Active)/2
	res <- fit$res
	
	while(any(G.old != G) & its < maxit){	
		d <- Inf
		while( d > tol & its < maxit){
			fit <- .C("grpsel", b0=as.double(b0), bigX = as.double(bigX),ssx=as.double(1), lambda= as.double(lambda), res=as.double(res),  
				Active=as.integer(t(Active)), n= as.integer(n), p=as.integer(p), k=as.integer(k), q=as.integer(q))
			d <- max(b0 - fit$b0)
			b0 <- fit$b0
			its <- its + 1
			
			res <- fit$res
			if(verbose) cat(".")
			}

		#Check to see if Active Set Changes:
		G.old <- G
		fit <- .C("grpsel", b0=as.double(b0), bigX = as.double(bigX), ssx=as.double(1), lambda= as.double(lambda), res=as.double(res),  
			Active=as.integer(t(all.Active)), n= as.integer(n), p=as.integer(p), k=as.integer(k), q=as.integer(q.max))
		
		b0 <- fit$b0
		b0.mat <-  matrix(b0 ,nrow=p,ncol=k*p,byrow=TRUE)
		G  <- b0.mat[,1:p] != 0
		Active <- AtoL(G)
		q <- length(Active)/2
		res <- fit$res
		
		if(verbose){
			ssres <- sum(fit$res^2)
			b2 <- b0.mat^2
			pen <- matrix(0,p,p)
			for(ipen in p*(1:k-1)) pen <- pen + b2[,ipen+1:p] + t(b2[,ipen+1:p])
			pen <- sum(sqrt(pen))/2
			cat("\n obj = ", lambda[1]*pen*(n-1)+(0.5)*ssres)
		}	
	}
	if(verbose) cat("\n")
	return(list(G=G,betas = matrix(b0 ,nrow=p,ncol=k*p,byrow=TRUE)))
}

.CovEM.Rcpp <- function(x, theta0, G, d, miss.group.unique, miss.group.counts, 
	miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, tol, maxiter, print.step){
	p <- ncol(x)
	res <- tryCatch( .Call("CovEM_Rcpp", x, nrow(x), p, theta0, G, d, 
		miss.group.unique, miss.group.counts, miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, tol, maxiter),
		"std::range_error" = function(e){
		conditionMessage( e ) } )
	iter <- res[1,1]
	eps <- res[2,2]
	mu <- res[3,]
	S <- res[1:p + 3,]
	## Reached max iteration messages
	#if( iter <= maxiter & print.step == 1) cat(paste("Reached convergence in", iter, "iterations.\n") )
	if(iter > maxiter) warning("Reached maximum number of iteration. No convergence is achieved.")
	list(mu=mu, S=S, iter=iter, eps=eps)
}


.CovEM.setparam <- function(p, mu0, S0)
{
	S0.uptri <- upper.tri(S0, diag=TRUE)
	theta <- c(-1,mu0, S0[S0.uptri])
	G.ind <- matrix(NA,p+1,p+1)
	G.ind[1,1] <- 1
	G.ind[1,-1] <- G.ind[-1,1] <- 1:p + 1
	G.ind[-1,-1][S0.uptri] <- 1:sum(S0.uptri) + p + 1
	G.ind[ lower.tri(G.ind)] <- t(G.ind)[ lower.tri(G.ind)]
	list( theta=theta, G.ind=G.ind)	
}
	

CovEM <- function(x, tol=0.001, maxiter=1000)
{
	xcall <- match.call()
	#if( !is.numeric(print.step) | print.step < 0 | print.step > 1) stop("argument 'print.step' must be: 0, 1.")

 	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")

	## drop all rows with missing values (!!) :
	n <- nrow(x); p <- ncol(x) 
	if( p < 2 ) stop("Column dimension of 'x' must be at least 2.")
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")
		
	## Rows with at least one observed
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	#if( length(not.ok) > 0 & print.step > 0) cat("Observations (rows): ", paste(not.ok, collapse=", "), 
	#	"\nare completely missing and will be dropped out from the estimation.\n")
	x.orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]	
	
	## get missing pattern and sorted matrix
	x_sort <- .sort.missing(x, x_nonmiss)

	## Check dimension 
	n <- nrow(x); p <- ncol(x)
	if(n <= p + 1)
		stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
	if(n < 2 * p)
		## p+1 < n < 2p
		warning("n < 2 * p, i.e., possibly too small sample size")

	## get initial estimate
	mu0 <- colMeans(x_sort$x, na.rm=T)
	S0 <- diag(apply(x_sort$x, 2, var, na.rm=T))
	x_sort <- c(x_sort, .CovEM.setparam(p, mu0, S0))

	## EM iteration
	res <- with(x_sort, .CovEM.Rcpp(x, theta, G.ind-1, length(theta), miss.group.unique, miss.group.counts, 
		miss.group.obs.col, miss.group.mis.col, miss.group.p, miss.group.n, tol, maxiter, print.step=0))
	
	## Compute pmd
	pmd.tmp <- .partial.mahalanobis.Rcpp( sweep(x_sort$x, 2, res$mu, "-"), res$S, x_sort$miss.group.unique, x_sort$miss.group.counts)
	pmd.tmp <- pmd.tmp[ x_sort$id.ro]
	pmd <- rep(NA, nrow(x.orig))
	pmd[ok] <- pmd.tmp	
	pmd.adj <- qchisq( pchisq( pmd, df=pp, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pp == p)] <- pmd[ which(pp==p) ] 
	
	## output results
	res <- new("CovRobMiss",
		call = xcall,
		S = res$S,
		mu = res$mu,
		estimator = "Maximum likelihood estimator", 
		x = x.orig,
		pmd = pmd,
		pmd.adj = pmd.adj,
		p = p,
		pu = pp )
	res
}





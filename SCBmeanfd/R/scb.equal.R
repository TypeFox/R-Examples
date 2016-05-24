scb.equal <-
function(x, y, bandwidth, level = .05, 
	scbtype = c("normal", "bootstrap", "both","no"), gridsize, keep.y = TRUE,
	nrep = NULL, nboot = NULL, parallel = c("no", "multicore", "snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
{
	caLL <- match.call()
	degree <- 1
	if (is.list(x)) {
		x1 <- x[[1]] ; x2 <- x[[2]]
	} else {
		x <- as.matrix(x)	
		if (ncol(x) == 1) {
			x1 = x2 = as.numeric(x) 
		} else if (ncol(x) == 2) {			
			x1 <- x[,1] ; x2 <- x[,2]
		} else stop("x has incorrect dimensions") 
	}	
	if(!isTRUE(all.equal(range(x1), range(x2)) )) 
		stop("x[[1]] and x[[2]] must have the same range")	
	spacing1 <- diff(x1)
	if (any(spacing1 < 0) || !isTRUE(all.equal(min(spacing1),max(spacing1))))
		stop("x[[1]] must be a uniform grid") 
	spacing2 <- diff(x2)
	if (any(spacing2 < 0) || !isTRUE(all.equal(min(spacing2),max(spacing2))))
		stop("x[[2]] must be a uniform grid") 
		
	y1 <- y[[1]]
	y2 <- y[[2]]
	n1 <- nrow(y1)
	n2 <- nrow(y2)		
	if (n1 < 2)
		stop("nrow(y[[1]]) must be greater than 1")
	if (n2 < 2)
		stop("nrow(y[[2]]) must be greater than 1")
	if (length(x1) != ncol(y1)) 
		stop("length(x[[1]]) and ncol(y[[1]]) do not match")
	if (length(x2) != ncol(y2)) 
		stop("length(x[[2]]) and ncol(y[[2]]) do not match")
	if (missing(gridsize))
		gridsize <- min(ncol(y1), ncol(y2))

	y1.hat <- apply(y1, 1, function(z) locpoly(x1, z, degree = degree, 
			  		bandwidth = bandwidth[1], gridsize = gridsize)$y)
	y2.hat <- apply(y2, 1, function(z) locpoly(x2, z, degree = degree, 
			 		bandwidth = bandwidth[2], gridsize = gridsize)$y)
	mu1.hat <- rowMeans(y1.hat)
	mu2.hat <- rowMeans(y2.hat)
	R1.hat <- cov(t(y1.hat)) / n1
	R2.hat <- cov(t(y2.hat)) / n2
	R.hat <- R1.hat + R2.hat
	std.error <- sqrt(diag(R.hat))
	test.stat <- max(abs(mu1.hat - mu2.hat) / std.error)

	p.norm = p.boot = NULL
	q.norm = q.boot = NULL
	lb.norm = ub.norm = lb.boot = ub.boot = NULL
	scbtype <- match.arg(scbtype)	

	if (scbtype %in% c("normal","both")) {
		eigcoR 	<- eigen(cov2cor(R.hat), TRUE)
		ncomp 	<- which(cumsum(eigcoR$values) > .99 * sum(eigcoR$values))[1]
		if(is.null(nrep)) nrep <- 2e4
		vars 	<- matrix(rnorm(ncomp * nrep), ncomp, nrep)
		M <- eigcoR$vectors[,1:ncomp] * rep(sqrt(eigcoR$values[1:ncomp]), each = gridsize)
		supnorm <- apply(abs(M %*% vars), 2, max)
		p.norm	<- 1 - ecdf(supnorm)(test.stat)
		q.norm	<- as.numeric(quantile(supnorm, 1-level)) 
		lb.norm	<- mu1.hat - q.norm * std.error
		ub.norm	<- mu1.hat + q.norm * std.error	
	}

	if (scbtype %in% c("bootstrap","both")) {
		if(is.null(nboot)) nboot <- 1e4		
		r <- t(cbind(y1.hat - mu1.hat, y2.hat - mu2.hat))
		ix1 <- 1:n1
		ix2 <- (n1+1):(n1+n2)
		boot.stat <- function(mat, ix) {	
			mat1 <- mat[ix[ix1],]
			mat2 <- mat[ix[ix2],]
			diffmu12.boot <- colMeans(mat1) - colMeans(mat2)
			sig12.boot <- sqrt( apply(mat1, 2, var) / n1 + apply(mat2, 2, var) / n2 )
			max(abs(diffmu12.boot) / sig12.boot)
			}
		supnorm <- boot(r, boot.stat, nboot, strata = rep(1:2, c(n1, n2)), 
			parallel = parallel, ncpus = ncpus, cl = cl)$t
		p.boot 	<- 1 - ecdf(supnorm)(test.stat)
		q.boot 	<- as.numeric(quantile(supnorm, 1 - level))
		lb.boot <- mu1.hat - q.boot * std.error
		ub.boot <- mu1.hat + q.boot * std.error	
	}

	result <- list( x = list(x1, x2), y = if(keep.y) y else NULL, call = caLL, model = NULL, par = NULL, nonpar = cbind(mu1.hat, mu2.hat), bandwidth = bandwidth, degree = degree, level = level, scbtype = scbtype, teststat = test.stat, pnorm = p.norm, pboot = p.boot, qnorm = q.norm, qboot = q.boot, normscb = cbind(lb.norm, ub.norm), bootscb = cbind(lb.boot, ub.boot), gridsize = gridsize, nrep = nrep, nboot = nboot )

	class(result) <- "SCBand"
	return(result)	
		
}

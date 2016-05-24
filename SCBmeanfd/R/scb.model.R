scb.model <-
function(x, y, model, bandwidth, level = .05, 
	scbtype = c("normal", "bootstrap", "both","no"), gridsize, keep.y = TRUE,
	nrep = NULL, nboot = NULL, parallel = c("no", "multicore", "snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
{
	caLL <- match.call()
	degree <- 1
	spacing <- diff(x)
	if (any(spacing < 0))
		stop("grid 'x' is not increasing") 
	if (!isTRUE(all.equal(min(spacing), max(spacing)))) 
		stop("grid 'x' is not uniform")		
	if (nrow(y) < 2)
		stop("nrow(y) > 1 is not TRUE")
	if (length(x) != ncol(y)) 
		stop("length(x) == ncol(y) is not TRUE")
	if (missing(gridsize))
		gridsize <- ncol(y)
		model <- as.matrix(model)
	if (nrow(model) != 1 && nrow(model) != length(x)) 
 		stop("'model' must be an integer or a matrix such that nrow(model) == length(x)")

	n <- nrow(y)
	if (length(model) == 1L && model == 0) {
		par.res <- t(y)
		par.mu.hat <- rep(0, gridsize)
	} 
	else if (length(model) == 1L && model > 0) {
		par.y.hat <- lm( t(y) ~ poly(x, degree = model))
		par.res <- residuals(par.y.hat)
		par.mu.hat <- rowMeans(fitted(par.y.hat))
	} 
	else {
		par.y.hat <- lm( t(y) ~ model - 1) 
		par.res <- residuals(par.y.hat)		    
		par.mu.hat <- rowMeans(fitted(par.y.hat))	
	}
	smooth.par.mu.hat <- locpoly(x, par.mu.hat, degree = degree, 
						   	bandwidth = bandwidth, gridsize = gridsize)$y
	nonpar.mu.hat <- locpoly(x, colMeans(y), degree = degree, 
							bandwidth = bandwidth, gridsize = gridsize)$y
	smooth.res <- apply(par.res, 2, function(z) locpoly(x, z, degree = 
					degree, bandwidth = bandwidth, gridsize = gridsize)$y)
	r <- rowMeans(smooth.res)
	sigma.hat <- apply(smooth.res, 1, sd)
	# scaled.smooth.res <- (smooth.res - r) / sigma.hat
	std.error <- sigma.hat / sqrt(n)
	test.stat <- max(abs(r / std.error))
	p.norm = p.boot = NULL
	q.norm = q.boot = NULL
	lb.norm = ub.norm = lb.boot = ub.boot = NULL
	scbtype <- match.arg(scbtype)	

	if (scbtype %in% c("normal","both")) {
		svd.sr <- svd((smooth.res - r)/ sigma.hat / sqrt(n-1), nv = 0)
		ncomp <- which(cumsum(svd.sr$d) > .99 * sum(svd.sr$d))[1]
		if(is.null(nrep)) nrep <- 2e4
		vars <- matrix(rnorm(ncomp * nrep), ncomp, nrep)
		M <- svd.sr$u[,1:ncomp] * rep(svd.sr$d[1:ncomp], each = gridsize)
		supnorm <- apply(abs(M %*% vars), 2, max)
		p.norm <- 1 - ecdf(supnorm)(test.stat)
		q.norm <- as.numeric(quantile(supnorm,1-level)) 
		lb.norm <- nonpar.mu.hat - q.norm * std.error
		ub.norm <- nonpar.mu.hat + q.norm * std.error	
	}

	if (scbtype %in% c("bootstrap","both")) {
		if(is.null(nboot)) nboot <- 5e3		
		boot.stat <- function(mat,ix) max(abs(colMeans(mat[ix,])) / apply(mat[ix,], 2, sd))
		supnorm <- boot(t(smooth.res - r), boot.stat, nboot, parallel = parallel, ncpus = ncpus, 
					cl = cl)$t * sqrt(n)
		p.boot <- 1 - ecdf(supnorm)(test.stat)
		q.boot <- as.numeric(quantile(supnorm,1-level))
		lb.boot <- nonpar.mu.hat - q.boot * std.error
		ub.boot <- nonpar.mu.hat + q.boot * std.error	
	}

	result <- list( x = x, y = if(keep.y) y else NULL, call = caLL, model = model, 
					par = smooth.par.mu.hat, nonpar = nonpar.mu.hat, 
					bandwidth = bandwidth, degree = degree, level = level, 
					scbtype = scbtype, teststat = test.stat,
					pnorm = p.norm, pboot = p.boot, 
					qnorm = q.norm, qboot = q.boot, 
					normscb = cbind(lb.norm, ub.norm), 
					bootscb = cbind(lb.boot, ub.boot), 
					gridsize = gridsize, nrep = nrep, nboot = nboot )

	class(result) <- "SCBand"
	return(result)	
		
}

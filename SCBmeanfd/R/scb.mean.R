scb.mean <-
function(x, y, bandwidth, level = .95, scbtype = c("normal", 
	"bootstrap", "both", "no"), gridsize, keep.y = TRUE, nrep = NULL, 
	nboot = NULL, parallel = c("no", "multicore", "snow"), ncpus = 
	getOption("boot.ncpus",1L), cl = NULL)
{
	caLL <- match.call()
	degree <- 1
	spacing <- diff(x)
	if (any(spacing < 0))
		stop("grid 'x' is not increasing") 
	if (!isTRUE(all.equal(min(spacing),max(spacing)))) 
		stop("grid 'x' is not uniform")		
	if (ncol(y) < 2)
		stop("ncol(y) > 1 is not TRUE")
	if (length(x) != ncol(y)) 
		stop("length(x) == ncol(y) is not TRUE")
	if (missing(gridsize))
		gridsize <- ncol(y)
	n <- nrow(y)
	N <- ncol(y)	
	y.hat <- apply(y, 1, function(z) locpoly(x, z, degree = degree, 
			 bandwidth = bandwidth, gridsize = gridsize)$y)
	mu.hat <- rowMeans(y.hat)
	sigma.hat <- apply(y.hat, 1, sd)
	std.error <- sigma.hat / sqrt(n)
	r <- y.hat - mu.hat  	
	lb.norm = ub.norm = lb.boot = ub.boot = q.norm = q.boot = NULL
	scbtype <- match.arg(scbtype)	

	if (scbtype %in% c("normal","both")) {
		svd.r <- svd(r / sigma.hat / sqrt(n-1), nv = 0)
		ncomp <- which(cumsum(svd.r$d) > .99 * sum(svd.r$d))[1]
		if(is.null(nrep)) nrep <- 2e4
		vars <- matrix(rnorm(ncomp * nrep), ncomp, nrep)
		M <- svd.r$u[,1:ncomp] * rep(svd.r$d[1:ncomp], each = gridsize)
		supnorm <- apply(abs(M %*% vars), 2, max)
		q.norm <- as.numeric(quantile(supnorm,level)) 
		lb.norm <- mu.hat - q.norm * std.error
		ub.norm <- mu.hat + q.norm * std.error	
	}
	if (scbtype %in% c("bootstrap","both")) {
		if(is.null(nboot)) nboot <- 5e3		
		boot.stat <- function(mat,ix) max(abs(colMeans(mat[ix,])/apply(mat[ix,], 2, sd)))
		supnorm   <- boot(t(r), boot.stat, nboot, parallel = parallel, ncpus = ncpus, cl = cl)$t
		q.boot 	  <- as.numeric(quantile(supnorm,level)) * sqrt(n)
		lb.boot   <- mu.hat - q.boot * std.error
		ub.boot   <- mu.hat + q.boot * std.error	
	}

	result <- list( x = x, y = if(keep.y) y else NULL, call = caLL, model = NULL,
					par = NULL, nonpar = mu.hat, bandwidth = bandwidth,
					degree = degree, level = level, scbtype = scbtype,
					teststat = NULL, pnorm = NULL, pboot = NULL, 
					qnorm = q.norm, qboot = q.boot,
					normscb = cbind(lb.norm, ub.norm), 
					bootscb = cbind(lb.boot, ub.boot), 
					gridsize = gridsize, nrep = nrep, nboot = nboot )
	class(result) <- "SCBand"
	return(result)
}

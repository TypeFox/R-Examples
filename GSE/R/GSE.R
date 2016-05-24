###########################################################
## Generalized S-Estimator
###########################################################
GSE <- function(x, tol=1e-4, maxiter=500, init=c("emve","qc","huber","imputed"), mu0, S0, ...)
{
	xcall <- match.call()

	## argument checks
	#method <- match.arg(method)
	method <- "bisquare"
	init <- match.arg(init)

	## check dat
	if(is.data.frame(x) | is.matrix(x))
		x <- data.matrix(x)
	else stop("Data matrix must be of class matrix or data.frame")
	
	p <- ncol(x)
	if( p >200 | p < 2 ) stop("Column dimension of 'x' must be in between 2 and 200.")
	
	## drop all rows with missing values (!!) :
	x_nonmiss <- is.na(x)*-1 + 1
	pp <- rowSums(x_nonmiss)
	pp_col <- colSums(x_nonmiss)
	
	## Cannot contain all obs with completely missing rows!!
	if( all(pp == 0) ) stop("All observations have missing values!")
	if( any(pp_col == 0) )stop("Data matrix cannot contain column(s) with completely missing data!")	
	
	ok <- which(pp > 0); not.ok <- which(pp == 0)
	x_orig <- x
	x <- x[ ok,]
	x_nonmiss <- x_nonmiss[ ok,]

	## reorder the data based on missingness
	x_sort <- .sort.missing(x, x_nonmiss)	

	## get initial estimate for all the EM calculation in EMVE subsampling
	EM.mu0 <- colMeans(x_sort$x, na.rm=T)
	EM.S0 <- diag(apply(x_sort$x, 2, var, na.rm=T))
	x_sort <- c(x_sort, .CovEM.setparam(p, EM.mu0, EM.S0))
	
	## dimension
	n <- nrow(x_sort$x); p <- ncol(x_sort$x)
	if(n <= p + 1)
		stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
	if(n < 2 * p)
		## p+1 < n < 2p
		warning("n < 2 * p, i.e., possibly too small sample size")

	if( xor(missing(mu0), missing(S0)) ) warning("Both 'mu0' and 'S0' must be provided. Default 'init' is used...")
	if( missing(mu0) || missing(S0) ){
		init.res <- switch( init,
			emve= { 
				with(x_sort, .emve.init(x, x_nonmiss, pu, n, p, theta, G.ind-1, length(theta),x.miss.group.match, 
					miss.group.unique, miss.group.counts, miss.group.obs.col, miss.group.mis.col, 
					miss.group.p, miss.group.n, ...))
				},
			qc ={res <- HuberPairwise(x, psi="sign", computePmd = FALSE); list(mu=res@mu, S=res@S) },
			huber = {res <- HuberPairwise(x, psi="huber", computePmd = FALSE, ...); list(mu=res@mu, S=res@S)},
			imputed = {ximp_simp <- .impute.simple(x, apply(x, 2, median, na.rm=TRUE)); 
					res <- rrcov::CovSest(ximp_simp, method=method);
					list(mu=res@center, S=res@cov) }
			)
		S0 <- init.res$S
		mu0 <- init.res$mu
	} 
	S0.chol <- tryCatch( chol(S0), error=function(e) NA)
	if( !is.matrix(S0.chol) )  stop("Estimated initial covariance matrix 'S0' is not positive definite.")

	## initiate GSE computation
	bdp <- 0.5
	print.step <- 0
	tol.scale <- 1e-9
	miter.scale <- 300
	res <- with(x_sort, .GSE.init(x, x_nonmiss, bdp, pu, n, p, miss.group.unique, miss.group.counts, mu0, S0, tol, 
		maxiter, tol.scale, miter.scale, print.step=print.step, method))
	
	## compute pmd
	pmd <- pmd.adj <- rep(NA, nrow(x_orig))
	pmd[ok] <- res$pmd[x_sort$id.ro]
	pmd.adj[ok] <- res$pmd.adj[x_sort$id.ro]	
	
	## new 2014-07-28
	wgts <- rep(NA, nrow(x_orig))
	wgts[ok] <- res$weights[x_sort$id.ro]
	wgtsp <- rep(NA, nrow(x_orig))
	wgtsp[ok] <- res$weightsprm[x_sort$id.ro]
	ximp <- matrix(NA, nrow(x_orig), ncol(x_orig))
	ximp[ok,] <- res$ximp[x_sort$id.ro,]
	
	res <- new("GSE",
		call = xcall,
		S = res$S,
		mu = res$mu,
		sc = res$stilde0,
		mu0 = mu0,
		S0 = S0, 
		iter = res$iter,
		eps = res$ep,
		estimator = "Generalized S-Estimator", 
		x = x_orig,
		ximp = ximp,
		weights = wgts,
		weightsp = wgtsp,
		pmd = pmd,
		pmd.adj = pmd.adj,
		p = p,
		pu = pp)
	res
}

## Assume the input data matrix is sorted using .sort.missing
.GSE.init <- function(x, x_nonmiss, bdp, pu, n, p, miss.group.unique, miss.group.counts,  mu0, S0, tol, maxiter, tol.scale, miter.scale, print.step, method)
{
	##########################################################################################
	## computing
	p.const.group <- rowSums(miss.group.unique)
	# if( method == "bisquare"){
		tuning.const.group <- .rho.bisquare.tune(p.const.group, bdp)
		res <- .GSE.Rcpp(x, matrix(mu0,1,p), S0, tol, maxiter, tol.scale, miter.scale, 
			miss.group.unique, miss.group.counts, tuning.const.group, print.step, bdp)
	# }
	# if( method == "rocke"){
		# tuning.const.group <- .rho.rocke.tune(p.const.group)
		# gamma.const.group <- pmin(1, qchisq(0.95, df=p.const.group)/p.const.group - 1)
		# res <- .GSE.rocke.Rcpp(x, matrix(mu0,1,p), S0, tol, maxiter, tol.scale, miter.scale, 
			# miss.group.unique, miss.group.counts, tuning.const.group, gamma.const.group, print.step, bdp)
	# }
	##########################################################################################
	## Include additional output other than mu and S output by .GSE.fixpt.init:
	## input data 'x', nonmissing variables 'pu', partial MD 'pmd', 
	## note: need to reorder x, pmd, pu
	
	## Compute pmd
	x.tmp <- sweep(x, 2, c(res$mu), "-")
	pmd <- .partial.mahalanobis.Rcpp(x.tmp, res$S, miss.group.unique, miss.group.counts)
	pmd.adj <- qchisq( pchisq( pmd, df=pu, log.p=T, lower.tail=F), df=p, log.p=T, lower.tail=F) 
	pmd.adj[ which( pu == p)] <- pmd[ which(pu==p) ]
	res$pmd <- pmd
	res$pmd.adj <- pmd.adj
	##########################################################################################
	## Output
	res
}






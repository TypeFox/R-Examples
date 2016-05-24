#
# gp class and method definitions
#

#' @keywords internal
#' @title gpsimple object, a simple Gaussian Process container with fields
#' @param x inputs
#' @param mean mean vector
#' @param cov covariance matrix
#' @param noisestd noise std vector
#' @param mll marginal log likelihood of data
#' @param x.obs observed input times
#' @param y.obs observed input values
#' @return a \code{gpsimple} object with same fields as parameters 
gpsimple = function(x, mean, cov, noisestd=NA, mll=NA, x.obs=NA,y.obs=NA) {
	gp = list('x'=x, 'mean'=mean, 'cov'=cov, 'noisestd'=noisestd, 'mll'=mll, x.obs=x.obs, y.obs=y.obs)
	class(gp) = 'gpsimple'
	return (gp)
}

#' @keywords internal
#' @title prints the gpsimple object
#' @param x the object
#' @param ... for compatibility
print.gpsimple = function(x, ...) {
	cat(length(x$x), 'timepoints from', min(x$x), 'to', max(x$x), '\n')
	cat('MLL:', x$mll, '\n')
}

#' @keywords internal
#' @title printss the gpsimple object
#' @description identical to \code{\link{print.gpsimple}}
#' @param object the object
#' @param ... for compatibility
summary.gpsimple = function(object, ...) {
	print(object)
}

sample.default = get('sample', mode='function')  # define existing sample() as default
sample = function(...) UseMethod('sample')       # make sample() generic
sample.gpsimple = function(obj, n=1) {             # define our own specialised sample()
	sample <- mvrnorm(n, obj$mean, obj$cov)
	return(sample)
}

sigmaupper = function(...) UseMethod('sigmaupper')
sigmaupper.gpsimple = function(obj, noisy=T, sigma=2) {
	res = obj$mean + sigma*sqrt(diag(obj$cov) + noisy*obj$noisestd^2)
	return (res)
}

sigmalower = function(...) UseMethod('sigmalower')
sigmalower.gpsimple = function(obj, noisy=T, sigma=2) {
	res = obj$mean - sigma*sqrt(diag(obj$cov) + noisy*obj$noisestd^2)
	return (res)
}

#' @export
#' @title prints the one-sample GP summary
#' @param x the estimated GP-object
#' @param ... for compatibility
print.gp = function(x, ...) {
	n = length(x$targets$x)
	cat('Gaussian process model for ', n, ' timepoints: (', sort(x$targets$x)[1], ', ', sort(x$targets$x)[2], ', ', (x$targets$x)[3], ', ..., ', sort(x$targets$x)[n-1], ', ', sort(x$targets$x)[n], ')\n\n', sep='')
	
	df = data.frame('MLL'=x$mll, 'EMLL'=x$emll, 'Avg posterior std'=mean(x$targets$pstd), 'Avg noise std'=mean(x$targets$noise), row.names='GP model')
	print(df, digits=3)
	
	cat('\nParameters:\n')
	cat(' sigma.f =', sprintf('%.2f', x$params$sigma.f), '\n')
	cat(' sigma.n =', sprintf('%.2f', x$params$sigma.n), '\n')
	cat('       l =', sprintf('%.2f', x$params$l), '\n')
	cat('    lmin =', sprintf('%.2f', x$params$lmin), '\n')
	cat('       c =', sprintf('%.2f', x$params$c), '\n')
}

#' @export
#' @title prints the one-sample GP summary
#' @description identical to \code{\link{print.gp}}
#' @param object the estimated GP-object
#' @param ... for compatibility
summary.gp = function(object, ...) {
	print(object)
}

# samples from the GP
sample.default = get('sample', mode='function')  # define existing sample() as default
sample = function(...) UseMethod('sample')       # make sample() generic
sample.gp = function(obj, N=1) {
	sample <- mvrnorm(N, obj$targets$pmean, obj$cov)
	return (sample)
}

sigmaupper = function(...) UseMethod('sigmaupper')
sigmaupper.gp = function(obj, noisy=T, sigma=2) {
	res = obj$targets$pmean + sigma*sqrt(obj$targets$pstd^2 + noisy*obj$targets$noisestd^2)
	return (res)
}

sigmalower = function(...) UseMethod('sigmalower')
sigmalower.gp = function(obj, noisy=T, sigma=2) {
	res = obj$targets$pmean - sigma*sqrt(obj$targets$pstd^2 + noisy*obj$targets$noisestd^2)
	return (res)
}


#' @export
#' @title prints the two-sample GP summary
#' @param x the estimated gppack-object containing ctrlmodel, casemodel and nullmodel
#' @param ... for compability
print.gppack = function(x, ...) {
	mlls  = c(x$ctrlmodel$mll, x$casemodel$mll, x$nullmodel$mll)
	pstds = c(mean(x$ctrlmodel$targets$pstd),mean(x$casemodel$targets$pstd),mean(x$nullmodel$targets$pstd))
	nstds = c(mean(x$ctrlmodel$targets$noise),mean(x$casemodel$targets$noise),mean(x$nullmodel$targets$noise))
	df = data.frame('MLL'=mlls, 'Avg posterior std'=pstds, 'Avg noise std'=nstds, row.names=c('Control model', 'Case model', 'Shared null model'))
	
	cat('Gaussian process models for case/control and shared null model\n\n')
	print(df, digits=3)
}

#' @export
#' @title prints the two-sample GP summary
#' @description identical to \code{\link{print.gppack}}
#' @param object the estimated gppack-object containing ctrlmodel, casemodel and nullmodel
#' @param ... for compability
summary.gppack = function(object, ...) {
	print(object)
}




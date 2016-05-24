#' Parametric Jackstraw
#'
#' Estimates statistical significance of association between variables and their latent variables, from a parametric jackstraw procedure.
#'
#' This function estimates statistical significance of association between variables and latent variables
#' using a parametric distribution of a noise term. A small number \code{s} of observed variables are replaced by
#' synthetic null variables generated from a specified distribution (such as Normal(0,1)).
#' After applying a latent variable estimation function on this newly generated matrix (with \code{s} synthetic nulls
#' and \code{m-s} intact observed variables), F-test statistics between estimated latent variables and \code{s} synthetic nulls
#' are called the jackstraw statistics. P-values are computed by comparing observed F-test statistics against \code{s*B} jackstraw statistics.
#'
#' Note that unless you have a strong reason to use a parametric distribution, it is advised to use the non-parametric jackstraw.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param FUN provide a function to estimate LVs. Must output \code{r} estimated LVs in a \code{n*r} matrix.
#' @param noise specify a parametric distribution to generate a noise term.
#' @param r a number of significant latent variables.
#' @param r1 a numeric vector of latent variables of interest.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param verbose a logical indicator as to whether to print the progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{jackstraw.parametric} returns a list consisting of
#' \item{p.value}{the \code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{the observed F-test statistics}
#' \item{null.stat}{the \code{s*B} null F-test statistics}
#'
#' @importFrom corpcor fast.svd
#' @import stats
#' @export jackstraw.parametric
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'
#' @seealso \link{jackstraw.FUN}
#' @seealso \link{jackstraw}
jackstraw.parametric = function(dat, FUN=function(x) fast.svd(x)$v[,1:r,drop=FALSE], noise=function(x) rnorm(x, mean=0, sd=1), r=NULL, r1=NULL, s=NULL, B=NULL, covariate=NULL, verbose=TRUE, seed=NULL) {
	if(!is.null(seed)) set.seed(seed)
  if(is.null(r)) stop("Must provide a number of latent variables, r.")
	m = dim(dat)[1]
	n = dim(dat)[2]
	if(is.null(s)) { s=round(m/10); message(paste0("A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=",s,".")); }
	if(is.null(B)) { B=round(m*10/s); message(paste0("A number of resampling iterations (B) is not specified: B=round(m*10/s)=",B,"."));}

	if(!is.null(FUN)) {
		FUN = match.fun(FUN)
		LV = FUN(dat)
		#if(is.null(r)) r = ncol(LV)
		if(r != ncol(LV)) stop(paste0("The number of latent variables ", r, "is not equal to the number of column(s) provided by ", FUN))
	} else {
		stop("Please provide a function to estimate latent variables.")
	}

	if(!is.null(noise)) {
		noise = match.fun(noise)
	} else {
		stop("Please provide a function to estimate latent variables.")
	}

	if(is.null(r1)) r1 = 1:r
	if(all(seq(r) %in% r1)) {
		# no adjustment LVs
		r0 = NULL
		ALV = NULL
		ALV.js = NULL
	} else {
		# r0 adjustment LVs
		r0 = seq(r)[-r1]
		ALV = LV[,r0, drop=FALSE]
		LV = LV[,r1, drop=FALSE]	## note that LV firstly contained the r latent variables; then reduced to the r1 latent variables of interest.
	}
	obs = FSTAT(dat=dat, LV=LV, ALV=ALV, covariate=covariate)$fstat

	# Estimate null association statistics
	null = matrix(0, nrow=s, ncol=B)

	if(verbose==TRUE) cat(paste0("\nComputating null statistics (", B," total iterations): "))
	for(i in 1:B){
		random.s = sample(1:m, size=s, replace=FALSE)
		s.nulls = matrix(noise(n*s), nrow=s, ncol=n)
		dat.js = dat
		dat.js[random.s,] = s.nulls

		LV.js = FUN(dat.js)
		if(!is.null(r0)) {
			ALV.js = LV.js[,r0, drop=FALSE]
			LV.js = LV.js[,r1, drop=FALSE]
		}
		null[,i] = FSTAT(dat=s.nulls, LV=LV.js, ALV=ALV.js, covariate=covariate)$fstat

		if(verbose==TRUE) cat(paste(i," "))
	}

  	p.value = cbind(getp(as.vector(obs), as.vector(null)))

	return(list(call=match.call(), p.value=p.value, obs.stat=obs, null.stat=null))
}

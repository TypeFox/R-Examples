#' Non-Parametric Jackstraw (Wrapper)
#'
#' Estimates statistical significance of association between variables and their latent variables (LVs).
#'
#' This is a wrapper for a few different functions using the jackstraw method.
#' Overall, it computes \code{m} p-values of association between \code{m} variables and their LVs.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of LVs from the observed data
#' and protects against an anti-conservative bias.
#'
#' For advanced use, one may consider computing association between variables and a subset of \code{r} estimated LVs.
#' For example, when there may be \code{r=3} significant PCs,
#' a user can carry out significance tests for the top two PCs (while adjusting for the third PC), by specifying \code{r1=c(1,2)} and \code{r=3}.
#'
#' Please take a careful look at your data and use appropriate graphical and statistical criteria
#' to determine a number of interesting/significant LVs, \code{r}.
#' It is assumed that \code{r} latent variables account for systematic variation in the data.
#'
#' For advanced usage, see \code{jackstraw.PCA}, \code{jackstraw.LFA}, and \code{jackstraw.FUN}.
#'
#' If \code{s} is not supplied, \code{s} is set to about 10\% of \code{m} variables.
#' If \code{B} is not supplied, \code{B} is set to \code{m*10/s}.
#'
#' @section Optional Arguments (see linked functions):
#' \describe{
#'   \item{s}{a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.}
#'   \item{B}{a number of resampling iterations.}
#'   \item{r1}{a numeric vector of latent variables (e.g., PCs) of interest. Not appropriate for all methods or functions.}
#'   \item{covariate}{a model matrix of covariates with \code{n} observations. Must include an intercept in the first column. Not appropriate for all methods and functions.}
#'   \item{verbose}{a logical specifying to print the computational progress. By default, \code{FALSE}.}
#'   \item{seed}{a seed for the random number generator.}
#' }
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param method a LV estimation method (by default, \code{"PCA"}). Use an optional argument \code{FUN} to specify a custom method.
#' @param FUN optionally, provide a specfic function to estimate LVs. Must output \code{r} estimated LVs in a \code{n*r} matrix.
#' @param r a number of significant latent variables.
#' @param ... optional arguemtns passed along to a specific jackstraw function.
#'
#' @return \code{jackstraw} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{\code{m} observed F-test statistics}
#' \item{null.stat}{\code{s*B} null F-test statistics}
#'
#' @importFrom corpcor fast.svd
#' @importFrom lfa lfa
#' @export jackstraw
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2013) Statistical significance of variables driving systematic variation in high-dimensional data Bioinformatics, 31(4): 545-554 \url{http://bioinformatics.oxfordjournals.org/content/31/4/545}
#'
#' @seealso \link{permutationPA} \link{jackstraw.PCA} \link{jackstraw.LFA} \link{jackstraw.FUN}
#'
#' @examples
#' set.seed(1234)
#' ## simulate data from a latent variable model: Y = BL + E
#' B = c(rep(1,50),rep(-1,50), rep(0,900))
#' L = rnorm(20)
#' E = matrix(rnorm(1000*20), nrow=1000)
#' dat = B %*% t(L) + E
#' dat = t(scale(t(dat), center=TRUE, scale=TRUE))
#'
#' ## apply the jackstraw
#' out = jackstraw(dat, r=1, method="PCA")
#'
#' ## Use optional arguments
#' ## For example, set s and B for a balance between speed of the algorithm and accuracy of p-values
#' \dontrun{
#' out = jackstraw(dat, r=1, s=10, B=1000, seed=5678)
#' }
jackstraw = function(dat, method="PCA", FUN=NULL, r=NULL, ...) {
	## for using Logistic Factors in linear regression (wrong..), use the following
	## out = jackstraw.FUN(dat, FUN = function(x) lfa(x, r)[,-(r+1),drop=FALSE], r = r, ...)
	message(paste("Using", method, "to estimate", r, "latent variables."))

	if(method == "PCA" || is.null(FUN)) {
		message("A linear regression with PCs, using a corpcor package.")
		out = jackstraw.PCA(dat, r = r, ...)
	} else if(method == "LFA") {
		message("A logistic regression with LFs, using a lfa package.")
		out = jackstraw.LFA(dat, FUN = function(x) lfa(x, r)[,,drop=FALSE], r = r, devR=FALSE, ...)
	} else if(method == "LFAcorpcor") {
		message("A logistic regression with LFs, using a corpcor package.")
		out = jackstraw.LFA(dat, FUN = function(x) lfa.corpcor(x, r)[,,drop=FALSE], r = r, devR=TRUE, ...)
	} else if(method == "custom" || !is.null(FUN)) {
		out = jackstraw.FUN(dat, FUN = FUN, r = r, ...)
	} else {
		stop("Wrong Inputs")
	}
	return(out)
}

#' Non-Parametric Jackstraw for Principal Component Analysis (PCA)
#'
#' Estimates statistical significance of association between variables and their principal components (PCs).
#'
#' This function computes \code{m} p-values of linear association between \code{m} variables and their PCs.
#' Its resampling strategy accounts for the over-fitting characteristics due to direct computation of PCs from the observed data
#' and protects against an anti-conservative bias.
#'
#' Provide the data matrix, with \code{m} variables as rows and \code{n} observations as columns.
#' Given that there are \code{r} significant PCs, this function tests for linear association between
#' \code{m} varibles and their \code{r} PCs.
#'
#' You could specify a subset of significant PCs that you are interested in (\code{PC}). If \code{PC} is given,
#' then this function computes statistical significance of association between \code{m} variables and \code{PC},
#' while adjusting for other PCs (i.e., significant PCs that are not your interest).
#' For example, if you want to identify variables associated with 1st and 2nd PCs,
#' when your data contains three significant PCs, set \code{r=3} and \code{PC=c(1,2)}.
#'
#' Please take a careful look at your data and use appropriate graphical and statistical criteria
#' to determine a number of significant PCs, \code{r}. The number of significant PCs depends on the data structure and the context.
#' In a case when you fail to specify \code{r}, it will be estimated from a permutation test (Buja and Eyuboglu, 1992)
#' using a function \link{permutationPA}.
#'
#' If \code{s} is not supplied, \code{s} is set to about 10\% of \code{m} variables.
#' If \code{B} is not supplied, \code{B} is set to \code{m*10/s}.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param r1 a numeric vector of principal components of interest. Choose a subset of \code{r} significant PCs to be used.
#' @param r a number (a positive integer) of significant principal components. See \link{permutationPA} and other methods.
#' @param s a number (a positive integer) of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number (a positive integer) of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param compute.obs a logical specifying to return observed statistics. By default, \code{TRUE}.
#' @param compute.null a logical specifying to return null statistics obtained by the jackstraw method. By default, \code{TRUE}.
#' @param compute.p a logical specifying to return p-values. By default, \code{TRUE}.
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a numeric seed for the random number generator.
#'
#' @return \code{jackstraw} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{\code{m} observed F-test statistics}
#' \item{null.stat}{\code{s*B} null F-test statistics}
#'
#' @importFrom corpcor fast.svd
#' @export jackstraw.PCA
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2013) Statistical Significance of Variables Driving Systematic Variation in High-Dimensional Data. arXiv:1308.6013 [stat.ME] \url{http://arxiv.org/abs/1308.6013}
#'
#' @seealso \link{jackstraw} \link{jackstraw.FUN} \link{permutationPA}
#'
#' @examples
#' set.seed(1234)
#' ## simulate data from a latent variable model: Y = BL + E
#' B = c(rep(1,50),rep(-1,50), rep(0,900))
#' L = rnorm(20)
#' E = matrix(rnorm(1000*20), nrow=1000)
#' dat = B %*% t(L) + E
#' dat = t(scale(t(dat), center=TRUE, scale=TRUE))
#'
#' ## apply the jackstraw
#' out = jackstraw.PCA(dat, r=1)
#'
#' ## Use optional arguments
#' ## For example, set s and B for a balance between speed of the algorithm and accuracy of p-values
#' \dontrun{
#' ## out = jackstraw.PCA(dat, r=1, s=10, B=1000, seed=5678)
#' }
jackstraw.PCA = function(dat, r1=NULL, r=NULL, s=NULL, B=NULL, covariate=NULL, compute.obs=TRUE, compute.null=TRUE, compute.p=TRUE, verbose=TRUE, seed=NULL) {
	m = dim(dat)[1]
	n = dim(dat)[2]
	if(is.null(s)) { s=round(m/10); message(paste0("A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=",s,".")); }
	if(is.null(B)) { B=round(m*10/s); message(paste0("A number of resampling iterations (B) is not specified: B=round(m*10/s)=",B,"."));}
	if(!is.null(seed)) set.seed(seed)
	if(is.null(r)) {
		warning("The number of significant PCs (r) is missing; this is strongly advised to determine r using appropriate statistical and graphical criteria.")
		r = permutationPA(dat=dat, threshold=.05, verbose=verbose)$r
		message(paste0("Permutation Parallel Analysis, with a threshold of 0.05, estimated r = ", r, "."))
	}
	if(!(r > 0 && r < n)) { stop("r is not in valid range between 1 and n-1."); }
	if(is.null(r1)) r1 = 1:r
	if(all(seq(r) %in% r1)) {
		# no adjustment LVs
		r0 = NULL
		ALV = NULL
	} else {
		r0 = seq(r)[-r1]
	}

	svd.dat = fast.svd(dat)
	LV = svd.dat$v[,r1, drop=FALSE]
	if(!is.null(r0)) ALV = svd.dat$v[,r0, drop=FALSE]

	# Calculate observed association statistics
	if(compute.obs==TRUE) {
	  obs = FSTAT(dat=dat, LV=LV, ALV=ALV, covariate=covariate)$fstat
	} else{
	  obs = NULL
	}

	# Estimate null association statistics
	if(compute.null==TRUE) {
  	null = matrix(0, nrow=s, ncol=B)
    ALV.js = NULL

  	if(verbose==TRUE) cat(paste0("\nComputating null statistics (", B," total iterations): "))
  	for(i in 1:B){
  		random.s = sample(1:m, size=s, replace=FALSE)
  		s.nulls = t(apply(dat[random.s, , drop=FALSE], 1, function(x) sample(x)))
  		dat.js = dat
  		dat.js[random.s,] = s.nulls

  		svd.dat.js = fast.svd(dat.js)
  		LV.js = svd.dat.js$v[,r1, drop=FALSE]
  		if(!is.null(r0)) ALV.js = svd.dat.js$v[,r0, drop=FALSE]
  		null[,i] = FSTAT(dat=s.nulls, LV=LV.js, ALV=ALV.js, covariate=covariate)$fstat

  		if(verbose==TRUE) cat(paste(i," "))
  	}
	} else {
	  null = NULL
	}

	if(compute.p==TRUE & compute.null==TRUE & compute.obs==TRUE) {
	  p.value = cbind(getp(as.vector(obs), as.vector(null)))
	} else {
	  p.value = NULL
	}

	return(list(call=match.call(), p.value=p.value, obs.stat=obs, null.stat=null))
}

#' Non-Parametric Jackstraw for Logistic Factor Analysis
#'
#' Estimates statistical significance of association between variables and their logistic factors (LFs).
#'
#' This function uses logistic factor analysis (LFA) from Wei et al. (2014). Particularly, dev in logistic regression
#' (the full model with \code{r} LFs vs. the intercept null model) is used to assess association.
#'
#' @param dat a genotype matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param FUN a function to use for LFA (by default, it uses the lfagen package)
#' @param devR use a R function to compute deviance. By default, FALSE (uses C++).
#' @param r a number of significant LFs.
#' @param r1 a numeric vector of LFs of interest (implying you are not interested in all \code{r} LFs).
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations. There will be a total of \code{s*B} null statistics.
#' @param covariate a data matrix of covariates with corresponding \code{n} observations (do not include an intercept term).
#' @param compute.obs a logical specifying to return observed statistics. By default, \code{TRUE}.
#' @param compute.null a logical specifying to return null statistics obtained by the jackstraw method. By default, \code{TRUE}.
#' @param compute.p a logical specifying to return p-values. By default, \code{TRUE}.
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{jackstraw} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their LFs}
#' \item{obs.stat}{\code{m} observed devs}
#' \item{null.stat}{\code{s*B} null devs}
#'
#' @importFrom corpcor fast.svd
#' @importFrom lfa lfa
#' @export jackstraw.LFA
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#'
#' @seealso \link{jackstraw} \link{jackstraw.FUN}
#'
#' @examples
#' set.seed(1234)
#' \dontrun{
#' ## simulate genotype data from a logistic factor model: drawing rbinom from logit(BL)
#' m=5000; n=100; pi0=.9
#' m0 = round(m*pi0)
#' m1 = m-round(m*pi0)
#' B = matrix(0, nrow=m, ncol=1)
#' B[1:m1,] = matrix(runif(m1*n, min=-.5, max=.5), nrow=m1, ncol=1)
#' L = matrix(rnorm(n), nrow=1, ncol=n)
#' BL = B %*% L
#' prob = exp(BL)/(1+exp(BL))
#'
#' dat = matrix(rbinom(m*n, 2, as.numeric(prob)), m, n)
#'
#' ## apply the jackstraw
#' out = jackstraw.LFA(dat, r=2)
#' }
jackstraw.LFA = function(dat, FUN = function(x) lfa(x, r)[,,drop=FALSE], devR=FALSE, r=NULL, r1=NULL, s=NULL, B=NULL, covariate=NULL, compute.obs=TRUE, compute.null=TRUE, compute.p=TRUE, verbose=TRUE, seed=NULL) {
	#FUN = function(x) lfa(x, r)[,,drop=FALSE]; r=1; s=NULL; B=NULL; covariate=NULL; verbose=TRUE; seed=NULL;
	if(!is.null(seed)) set.seed(seed)

	m = dim(dat)[1]
	n = dim(dat)[2]
	if(is.null(s)) { s=round(m/10); message(paste0("A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=",s,".")); }
	if(is.null(B)) { B=round(m*10/s); message(paste0("A number of resampling iterations (B) is not specified: B=round(m*10/s)=",B,"."));}

	if(!(r > 0 && r < n)) { stop("r is not in valid range between 1 and n-1."); }
	if(is.null(r1)) r1 = 1:r
	if(all(seq(r) %in% r1)) {
		# no adjustment LVs
		r0 = NULL
		LFr0 = NULL
	} else {
		r0 = seq(r)[-r1]
	}

	if(!is.null(FUN)) {
		FUN = match.fun(FUN)
	} else {
		stop("Please provide a function to estimate latent variables.")
	}

	## note that LFr has an intercept term as the last column
	## LFr1 and LFr0 (subsetting LFr) do not inherit/have an intercept term
	LFr = FUN(dat)
	LFr1 = LFr[,r1, drop=FALSE]
	if(r != ncol(LFr)) stop(paste0("The number of latent variables ", r, "is not equal to the number of column(s) provided by ", FUN))
	if(!is.null(r0)) LFr0 = LFr[,r0, drop=FALSE]

	if(compute.obs == TRUE) {
  	if(devR == FALSE) {
  		## uses a deviance computation function in lfagen
  		obs = devdiff(dat, LF_alt=cbind(LFr, covariate), LF_null=cbind(LFr0, matrix(1, n, 1), covariate))
  	} else {
  		## uses a deviance computation function from R base
  		obs = dev.R(dat, LFr1=cbind(LFr1, covariate), LFr0=cbind(LFr0, covariate))
  	}
	} else {
	  obs = NULL
	}

	# Estimate null association statistics
	if(compute.null == TRUE) {
  	null = matrix(0, nrow=s, ncol=B)
    LFr0.js = NULL

  	if(verbose==TRUE) cat(paste0("\nComputating null statistics (", B," total iterations): "))
  	for(i in 1:B){
  		random.s = sample(1:m, size=s, replace=FALSE)
  		s.nulls = t(apply(dat[random.s, , drop=FALSE], 1, function(x) sample(x)))
  		dat.js = dat
  		dat.js[random.s,] = s.nulls

  		LFr.js = FUN(dat.js)
  		LFr1.js = LFr.js[,r1, drop=FALSE]
  		if(!is.null(r0)) LFr0.js = LFr.js[,r0, drop=FALSE]

  		if(devR==FALSE) {
  			## uses a deviance computation function in lfagen
  			null[,i] = devdiff(s.nulls, LF_alt=cbind(LFr.js, covariate), LF_null=cbind(LFr0.js, matrix(1, n, 1), covariate))
  		} else {
  			## uses a deviance computation function from R base
  			null[,i] = dev.R(s.nulls, LFr1=cbind(LFr1.js, covariate), LFr0=cbind(LFr0.js, covariate))
  		}

  		if(verbose==TRUE) cat(paste(i," "))
  	}
	} else{
	  null = NULL
	}

  if(compute.p==TRUE & compute.null==TRUE & compute.obs==TRUE) {
    p.value = cbind(getp(as.vector(obs), as.vector(null)))
  } else {
    p.value = NULL
  }

	return(list(call=match.call(), p.value=p.value, obs.stat=obs, null.stat=null))
}

#' Non-Parametric Jackstraw for a Custom Function
#'
#' Estimates statistical significance of association between variables and their latent variables, estimated using a custom function.
#'
#' @param dat a data matrix with \code{m} rows as variables and \code{n} columns as observations.
#' @param FUN optionally, provide a specfic function to estimate LVs. Must output \code{r} estimated LVs in a \code{n*r} matrix.
#' @param r a number of significant latent variables.
#' @param r1 a numeric vector of latent variables of interest.
#' @param s a number of ``synthetic'' null variables. Out of \code{m} variables, \code{s} variables are independently permuted.
#' @param B a number of resampling iterations.
#' @param covariate a model matrix of covariates with \code{n} observations. Must include an intercept in the first column.
#' @param compute.obs a logical specifying to return observed statistics. By default, \code{TRUE}.
#' @param compute.null a logical specifying to return null statistics obtained by the jackstraw method. By default, \code{TRUE}.
#' @param compute.p a logical specifying to return p-values. By default, \code{TRUE}.
#' @param verbose a logical specifying to print the computational progress.
#' @param seed a seed for the random number generator.
#'
#' @return \code{jackstraw} returns a list consisting of
#' \item{p.value}{\code{m} p-values of association tests between variables and their principal components}
#' \item{obs.stat}{\code{m} observed statistics}
#' \item{null.stat}{\code{s*B} null statistics}
#'
#' @importFrom corpcor fast.svd
#' @export jackstraw.FUN
#' @author Neo Christopher Chung \email{nchchung@@gmail.com}
#' @references Chung and Storey (2013) Statistical Significance of Variables Driving Systematic Variation in High-Dimensional Data. arXiv:1308.6013 [stat.ME] \url{http://arxiv.org/abs/1308.6013}
#'
#' @seealso \link{jackstraw}
#'
#' @examples
#' set.seed(1234)
#' ## simulate data from a latent variable model: Y = BL + E
#' B = c(rep(1,50),rep(-1,50), rep(0,900))
#' L = rnorm(20)
#' E = matrix(rnorm(1000*20), nrow=1000)
#' dat = B %*% t(L) + E
#' dat = t(scale(t(dat), center=TRUE, scale=TRUE))
#'
#' ## apply the jackstraw with the svd as a function
#' out = jackstraw.FUN(dat, FUN = function(x) svd(x)$v[,1,drop=FALSE], r=1, s=100, B=50)
jackstraw.FUN = function(dat, FUN, r=NULL, r1=NULL, s=NULL, B=NULL, covariate=NULL, compute.obs=TRUE, compute.null=TRUE, compute.p=TRUE, verbose=TRUE, seed=NULL) {
	if(!is.null(seed)) set.seed(seed)

	m = dim(dat)[1]
	n = dim(dat)[2]
	if(is.null(s)) { s=round(m/10); message(paste0("A number of null variables (s) to be permuted is not specified: s=round(0.10*m)=",s,".")); }
	if(is.null(B)) { B=round(m*10/s); message(paste0("A number of resampling iterations (B) is not specified: B=round(m*10/s)=",B,"."));}

	if(!is.null(FUN)) {
		FUN = match.fun(FUN)
		LV = FUN(dat)
		if(is.null(r)) r = ncol(LV)
		if(r != ncol(LV)) stop(paste0("The number of latent variables ", r, "is not equal to the number of column(s) provided by ", FUN))
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

	if(compute.obs==TRUE) {
  	obs = FSTAT(dat=dat, LV=LV, ALV=ALV, covariate=covariate)$fstat
	} else {
	  obs = NULL
	}

	# Estimate null association statistics
	if(compute.null==TRUE) {
  	null = matrix(0, nrow=s, ncol=B)

  	if(verbose==TRUE) cat(paste0("\nComputating null statistics (", B," total iterations): "))
  	for(i in 1:B){
  		random.s = sample(1:m, size=s, replace=FALSE)
  		s.nulls = t(apply(dat[random.s, , drop=FALSE], 1, function(x) sample(x)))
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
	} else {
	  null = NULL
	}

	if(compute.p==TRUE & compute.null==TRUE & compute.obs==TRUE) {
	  p.value = cbind(getp(as.vector(obs), as.vector(null)))
	} else {
	  p.value = NULL
	}

	return(list(call=match.call(), p.value=p.value, obs.stat=obs, null.stat=null))
}


# R package: generalized meta-analysis (gmeta), file gmeta.r 
# R interface to meta-analysis based on combine confidence distributions (CDs)

# *****************************************************************************
#    main
# *****************************************************************************

# main function
gmeta <- function(gmi, gmi.type = c('pivot', 'cd', 'pvalue', '2x2'),
		method = c('fixed-mle',
		        'fixed-robust1', 'fixed-robust2', 'fixed-robust2(sqrt12)',
				'random-mm', 'random-reml', 'random-tau2',
				'random-robust1', 'random-robust2', 'random-robust2(sqrt12)',
				'fisher', 'normal', 'stouffer', 'min', 'tippett', 'max', 'sum',
				'MH', 'Mantel-Haenszel', 'Peto', 'exact1', 'exact2'),
		linkfunc = c('inverse-normal-cdf', 'inverse-laplace-cdf'), 
		weight = NULL, study.names = NULL, gmo.xgrid = NULL, ci.level = 0.95, 
		tau2 = NULL, mc.iteration = 10000, eta = 'Inf', verbose = FALSE, 
		report.error = FALSE)
{ 
	UseMethod('gmeta')
}

# main function
gmeta.default <- function(gmi, # generalized meta-analysis input: vector of p-values, data.frame/matrix of mean/sd, list of cd-functions or matrix of 2x2 tables.
		# parameters specifying input properties
		gmi.type = c('pivot',  # type of input (gmi), if model-based meta-analysis, input data.frame/matrix of mean/sd
		             'cd',     # type of input (gmi), if model-based meta-analysis, input list of cd-funcitons (pointers to function)
					 'pvalue', # type of input (gmi), if combine p-values, input vector of p-values
					 '2x2'),   # type of input (gmi), if combine 2x2 tables, input matrix of 2x2 tables in format c(rx-event_treatment_grounp, nx-number_observation_treatment_group, ry-event_control_group, ny-number_observation_treatment_group) for each row of the matrix
		method   = c('fixed-mle',      # method of combination, if model-based meta-analysis, fixed-effects model with mle
		             'fixed-robust1',  # method of combination, if model-based meta-analysis, fixed-effects model with robust method 1 - adaptive-weighting (xie2012section4)
					 'fixed-robust2',  # method of combination, if model-based meta-analysis, fixed-effects model with robust method 2 - m-estimating-equation (xie2012section4)
					 'fixed-robust2(sqrt12)', # method of combination, if model-based meta-analysis, fixed-effects model with robust method 2 - m-estimator-simplified-variance (xie2012section4)
				     'random-mm',      # method of combination, if model-based meta-analysis, random-effects model with moment estimator
					 'random-reml',    # method of combination, if model-based meta-analysis, random-effects model with restricted-mle estimator
					 'random-tau2',    # method of combination, if model-based meta-analysis, random-effects model with heterogeneity specified by user
				     'random-robust1', # method of combination, if model-based meta-analysis, random-effects model with robust method 1 - adaptive-weighting (xie2012section4)
					 'random-robust2', # method of combination, if model-based meta-analysis, random-effects model with robust method 2 - m-estimating-equation (xie2012section4)
					 'random-robust2(sqrt12)', # method of combination, if model-based meta-analysis, random-effects model with robust method 2 - m-estimator-simplified-variance (xie2012section4)
				     'fisher',   # method of combination, if combine p-value vector, see fisher1932, bahadur optimal
					 'normal',   # method of combination, if combine p-value vector, see stouffer(1949)
					 'stouffer', # method of combination, if combine p-value vector, see stouffer1949, alias 'normal'
					 'min',      # method of combination, if combine p-value vector, see tippett(1931)
					 'tippett',  # method of combination, if combine p-value vector, see tippett1931, alias 'min'
					 'max',      # method of combination, if combine p-value vector, see marden(1991)
					 'sum',      # method of combination, if combine p-value vector, see marden(1991)
				     'MH',   # method of combination, if cobmine 2x2 tables, see Mantel-Haenszel method
					 'Mantel-Haenszel', # method of combination, if cobmine 2x2 tables, see Mantel-Haenszel method, alias 'MH'
					 'Peto', # method of combination, if cobmine 2x2 tables, see Peto's method
					 'exact1',  # method of combination, if cobmine 2x2 tables, see liu2012exact
					 'exact2'), # method of combination, if cobmine 2x2 tables, see tian2009exact
		linkfunc = c('inverse-normal-cdf',   # 'link' function, F0^{-1}(\cdot), see xie2012section2 
		             'inverse-laplace-cdf'), # 'link' function, F0^{-1}(\cdot), see xie2012section2
		weight   = NULL, # study-specific weight, in default inverse-of-adjusted-variance
		study.names = NULL, # names of studies, in default, study1, study2, etc.
		# parameters specifying output properties
		#gmo.limit= NULL, # range of generalized meta-analysis output, in default, c(-1,1) with gmo.n.pts=2001, seq(from=-1,to=1,by=0.001)
		#gmo.n.pts= 2001, # number of points in the range of generalized meta-analysis output, suppressed if gmo.limit is given as a vector 
		gmo.xgrid= NULL,  # gridding positions where generalized meta-analysis output (fixed- or random-effects model-based) evaluating the combined and individual CDs. gmo.xgrids suppress gmo.limit and gmo.n.pts, in default, gmo.xgrid = seq(-1,1,by=0.001).
		ci.level = 0.95,  # confidence level
		# parameters specifying model-based meta-analysis (random-effects model) heterogeneity
		tau2     = NULL,  # numeric value or method used to estimate heterogeneity, in default, REML. (only when method=='random-tau2')
		# parameters specifying properties when use liu2012exact (exact1) method to combine 2x2 tables
		mc.iteration = 10000, # monte-carlo iteration (see liu2012exact)
		# parameters specifying properties when use tian2009exact (exact2) method to combine 2x2 tables
		eta      = 'Inf', # seq(0.05, 0.95, length=20), # number of confidence-level taken to make confidence interval for later combination (see tian2009exact)
		# parameters specifying whether report detailed inforamtion
		verbose  = FALSE, # whether detailed inforamtion reported during meta-analysis
		# parameters specifying whether report error of coverage probability when use liu2012exact/tian2009exact method to combine 2x2 tables
		report.error = FALSE) # minor error message and exact error of coverage probability of the computed confidence interval for 2x2-"exact1"&"exact2" method
{ 
	#check
	mf          <- match.call()
	gmi.type    <- match.arg(gmi.type)
	gmi         <- gmeta.check.gmi(gmi, gmi.type, verbose)
	method      <- match.arg(method)
	linkfunc    <- match.arg(linkfunc)
	weight      <- gmeta.check.weight(gmi, gmi.type, linkfunc, weight, verbose)
	gmo.xgrid   <- gmeta.check.gmo.xgrid(gmi, gmi.type, gmo.xgrid)
	study.names <- gmeta.check.study.names(gmi, study.names)
	#gmo.n.pts <- round(ifelse(is.numeric(gmo.n.pts)&&(gmo.n.pts>1), gmo.n.pts, 2001))
	#if(! ((!is.null(gmo.limit))&&is.numeric(gmo.limit)&&(length(gmo.limit)==2)) ) { gmo.limit = c(-1,1) }
	ci.level  <- ifelse(is.numeric(ci.level)&&(ci.level>=0)&&(ci.level<=1), ci.level, 0.95)
	mc.iteration <- round(ifelse(is.numeric(mc.iteration)&&(mc.iteration>1), mc.iteration, 10000))
	#if(!(is.numeric(eta)&&(min(eta)>=0)&&(max(eta)<=1))) { eta = seq(0.05, 0.95, length=20) } # allow eta non-numeric as K\to\Infty
	# meta-analysis
	if (gmi.type == 'pvalue') {
		gmo <- gmeta.p(gmi, method)
	} else if (gmi.type == 'cd' || gmi.type == 'pivot') {
		gmo <- gmeta.m(gmi, gmi.type, method, linkfunc, weight, gmo.xgrid, ci.level, tau2, verbose)
	} else if (gmi.type == '2x2') { # combine 2x2 tables
		gmo <- gmeta.e(gmi, method, weight, gmo.xgrid, ci.level, mc.iteration, eta, verbose, report.error)
	} else {
		stop("gmi.type must be 'cd', 'pivot', 'pvalue' or '2x2'.")
	}
	# post processing
	gmo$call   <- match.call()
	gmo$input  <- gmi
	gmo$alpha  <- 1 - ci.level
	gmo$study.names <- study.names
	# registr S3 class
	if ( gmi.type=='pvalue' ) {
		class(gmo) <- c('gmeta.p', 'gmeta')
	} else if ( gmi.type == 'cd' || gmi.type == 'pivot' ) {
		class(gmo) <- c('gmeta.m', 'gmeta')
	} else if (gmi.type == '2x2') {
		class(gmo) <- c('gmeta.e', 'gmeta')
	} else {
		stop("gmi.type must be 'cd', 'pivot', 'pvalue' or '2x2'.")
	}
	# return
	return(gmo)
}




# preprocessing

## preprocessing - checks
# *****************************************************************************
#    # checking - input, output, parameters formats, etc.
# *****************************************************************************
## main [gmeta.check.x()]

### check gmi - gmeta input formats
### update v2.0: 2x2 table order is changed - before x y M N, now x-M y-N.
gmeta.check.gmi <- function(gmi, gmi.type, verbose) {
	# check by gmi.type
	if ( gmi.type == 'cd' ) {
		# verbose
		if ( verbose ) {
			cat('\ninput CDs must composite of a list.\n') #commentOut
		}
		# check
		if ( !is.list(gmi) ) {
			stop('input CDs must composite of a list.')
		}
	} else if ( gmi.type == 'pivot' ) {
		# verbose
		if ( verbose ) {
			cat('\ninput pivots must be in form of of either a matrix or a data.frame of 2 columns.\n') #commentOut
		}
		# check
		if ( ! (is.matrix(gmi) || is.data.frame(gmi)) ) {
			stop('input pivots must be in form of either a matrix or a data.frame.')
		}
		if ( dim(gmi)[2] != 2 ) {
			stop('input pivots must be in form of of either a matrix or a data.frame of 2 columns.')
		}
	} else if ( gmi.type == 'pvalue' ) {
		# verbose
		if ( verbose ) {
			cat('\ninput p-values must be in form of a vector with all element in-between 0 and 1.\n') #commentOut
		}
		# check
		if ( !is.vector(gmi) ) {
			stop('input p-values must be in form of a vector.')
		}
		if ( any( gmi < 0 | gmi > 1 ) ) {
			stop('input p-values must be between 0 and 1.')
		}
	} else if ( gmi.type == '2x2' ) {
		# verbose
		if ( verbose ) {
			cat('\ninput 2x2 tables should be in format: r_group_1 n_group_1 r_group_2 n_group_2.\n') #commentOut
		}
		# check
		if ( ! (is.matrix(gmi) || is.data.frame(gmi)) ) {
			stop('input 2x2 tables must be in form of either a matrix or a data.frame.')
		}
		if ( dim(gmi)[2] != 4 ) {
			stop('input 2x2 tables must composites of either a matrix or a data.frame of 4 columns.')
		}
		if ( any( gmi[,1] > gmi[,2] | gmi[,3] > gmi[,4] ) ) {
			stop('input 2x2 tables should be in format: r_group_1 n_group_1 r_group_2 n_group_2 - second and fourth column should be larger than first and third column, respectively.')
		}
		if ( any( gmi < 0 ) ) {
			stop('input data should be nonnegative.')
		}
	} else {
		stop('gmi.type not recognize.')
	}
	# return
	return(gmi)
}

### check weight - weight assigned to each study
### update v2.0: weights: need to distinguish whether use weight=1/s^2 (gaussian) or weight=1(laplacian, DE, warning if not), keep NULL set later.
gmeta.check.weight <- function(gmi, gmi.type, linkfunc, weight, verbose) {
	if ( gmi.type == 'cd' || gmi.type == 'pivot' ) {
		# number of studies
		if ( gmi.type == 'cd' ) {
			nn1 = length(gmi)
		} else { # gmi.type == 'pivot'
			nn1 = dim(gmi)[1]
		}
		if ( (!is.null(weight)) && (is.vector(weight)) ) {
			nn2 = length(weight) # length of weight vector
			if (nn2 != nn1) {
				stop('weight must be assigned to each study.')
			}
			if (!match(typeof(weight), c('double', 'integer', 'logical'))) {
				stop('weight must be a vector comprised by numeric numbers.')
			}
		} else if ( (!is.null(weight)) && (!is.vector(weight)) ) {
			stop('weight must be a vector comprised by numeric numbers.')
		} else {
			# verbose
			if ( verbose ) {
				cat('\nweight is null - use default weight which depends on method.\n') #commentOut
			} else {
				#cat('\nweight is null - use default weight which depends on method.\n') #commentOut
			}
		}
	} # note: weight=NULL is kept.
	else if ( gmi.type == 'pvalue' ) { 
		if ( is.vector(gmi) ) {
			if ( is.null(weight) ) {
				# verbose
				if ( verbose ) {
					cat('\nweight is null - use default weight all one.\n') #commentOut
				}
				weight=rep(1, length(gmi))
			} else if ( is.numeric(weight) ) {
				if ( length(weight) != length(gmi) ) {
					stop('weight must be assigned to each trial.')
				}
			} else {
				stop('weight must be a vector comprised by numeric numbers.')
			}
		} else {
			stop('input for p-value combination must be a vector.')
		}
	}
	else if ( gmi.type == '2x2' ) {
		#if (is.null(weight)) {
			# weight = rep(1, dim(gmi)[1])
		#}
		if ( !is.null(weight) ) {
			stopifnot( is.numeric(weight), length(weight) == dim(gmi)[1] )
		} # note: weight=NULL is kept.
	}
	else {
		stop('input type must be cd, pivot, pvalue or 2x2.')
	}
	# return
	return(weight)
}

### checking gmo.xgrid: cd, pivot, exact
#####
gmeta.check.gmo.xgrid <- function(gmi, gmi.type, gmo.xgrid) {
	if (gmi.type == 'cd') {
		gmo.xgrid <- gmeta.check.gmo.xgrid.cd(gmi, gmo.xgrid)
	}
	else if (gmi.type == 'pivot') {
		gmo.xgrid <- gmeta.check.gmo.xgrid.pivot(gmi, gmo.xgrid)
	}
	else if (gmi.type == 'pvalue') { 
		gmo.xgrid <- gmeta.check.gmo.xgrid.pvalue(gmi, gmo.xgrid)
	}
	else if (gmi.type == '2x2') {
		gmo.xgrid <- gmeta.check.gmo.xgrid.2x2(gmi, gmo.xgrid)
	}
	else {
		stop('input type must be cd, pivot, pvalue or 2x2!')
	}
	return(gmo.xgrid)
}
#####
gmeta.check.gmo.xgrid.cd <- function(gmi, gmo.xgrid) {
	if( is.null(gmo.xgrid) ) {
		nn1 <- length(gmi)
		gmi.x <- list()
		for (i in seq(1, nn1)) {
			gmi.x[[i]] <- gmi[[i]][,1]
		}
		xl <- min(unlist(gmi.x))
		xu <- max(unlist(gmi.x))
		gmo.xgrid <- seq(xl, xu, by=0.001)
	} else {
		stopifnot( is.vector(gmo.xgrid), is.numeric(gmo.xgrid) )
	}
	return(gmo.xgrid)
}
#####
gmeta.check.gmo.xgrid.pivot <- function(gmi, gmo.xgrid) {
	if ( is.null(gmo.xgrid) ) {
		gmo.xgrid <- seq(from=-1, to=1, by=0.001)
	} else {
		stopifnot( is.vector(gmo.xgrid), is.numeric(gmo.xgrid) )
	}
	return(gmo.xgrid)
}
#####
gmeta.check.gmo.xgrid.pvalue <- function(gmi, gmo.xgrid) {
	if ( !is.null(gmo.xgrid) ) {
		warning('gmo.xgrid is not used in p-value combination.')
	}
	return(NULL)
}
#####
gmeta.check.gmo.xgrid.2x2 <- function(gmi, gmo.xgrid){
	if( is.null(gmo.xgrid) ) {
		gmo.xgrid = seq(from=-1, to=1, by=0.001)
	} else {
		stopifnot( is.vector(gmo.xgrid), is.numeric(gmo.xgrid), length(gmo.xgrid) > 2 )
	}
	return(gmo.xgrid)
}

### checking study.names
gmeta.check.study.names <- function(gmi, study.names) {
	# nn1 - number of study
	if ( is.matrix(gmi) || is.data.frame(gmi) ) {
		nn1 = dim(gmi)[1]
	} else if ( is.list(gmi) || is.vector(gmi) ) {
		nn1 = length(gmi)
	} else {
		stop('impossible')
	}
	# nn2 - number of study.names
	if ( is.null(study.names) ) {
		# study.names
		if ( is.list(gmi) && !is.data.frame(gmi) && !is.matrix(gmi) ) { # just a list
			study.names <- names(gmi)
		} else if ( is.data.frame(gmi) || is.matrix(gmi) ) { # R - is.list(data.frame)=TRUE
			study.names <- rownames(gmi)
		} else {
			study.names = paste('study-', formatC(c(1:nn1),width=ceiling(log10(nn1)),format='d',flag='0'), sep='')
		}
		# study.names - again
		if ( is.null(study.names) || (study.names == rep(1,nn1)) ) {
			study.names = paste('study-', formatC(c(1:nn1),width=ceiling(log10(nn1)),format='d',flag='0'), sep='')
		}
	} else if ( is.vector(study.names) ) {
		nn2 = length(study.names)
		if ( !(nn1 == nn2) ) {
			stop('study.names must assign each study a name.')
		}
	} else {
		stop('if not NULL, study.names must be a vector that assigns each study a name.')
	}
	# return
	return(study.names)
}

# preprocessing[done]




# A unifying framework for meta-analysis 

## meta-analysis - combine p-value vector
# *****************************************************************************
#    gmeta.p() - combination of p-value vector
# *****************************************************************************
## main [gmeta.p()]
gmeta.p <- function(gmi, method) {
	cmbd.pvalue <- Cpvaluecombine(gmi, method)
	gmeta.cmbd <- list(individual.pvalues=gmi, method=method, combined.pvalue=cmbd.pvalue)
	return(gmeta.cmbd)
}
### combine: p-value vector
### Cpvaluecombine: implemented in C
Cpvaluecombine <- function(RpVec, Rmethod) {
	stopifnot(is.numeric(RpVec), is.character(Rmethod));
	stopifnot(Rmethod %in% c('fisher', 'normal', 'stouffer', 'min', 'tippett', 'max', 'sum'));
	.Call('pvaluecombine', # pvaluecombine() in C
	RpVec,   # generatlized meta-analysis input: vector of p-values.
	Rmethod, # generatlized meta-analysis method: 'fisher', 'normal', 'stouffer', 'min', 'tippett', 'max', 'sum'. 
	package  = 'pvaluecombine')
}
### combine: p-value vector
### R.pvalue.combine: implemented in R (deprecated)
R.pvalue.combine <- function(gmi, method) {
	k = length(gmi)
	if (method == 'fisher') {
		cmbd.pvalue <- 1 - pchisq(-2*sum(log(gmi)), df=2*k)  
	}
	else if (method == 'normal') {
		cmbd.pvalue <- pnorm(sum(qnorm(gmi))/sqrt(k))
	}
	else if (method == 'stouffer') {
		# alias for normal
		cmbd.pvalue <- pnorm(sum(qnorm(gmi))/sqrt(k))
	}
	else if (method == 'tippett') {
		cmbd.pvalue <-  1-(1-min(gmi))^k
	}
	else if (method == 'max') {
		cmbd.pvalue <- max(gmi)^k
	}
	else if (method == 'min') {
		# alias for 'tippett'
		cmbd.pvalue <-  1-(1-min(gmi))^k
	}
	else if (method == 'sum') {
		if (k <= 30) {
			cmbd.pvalue <- psumunif(sum(gmi), k)
		}
		else {
			cmbd.pvalue <- pnorm(sum(gmi), mean = k/2, sd = sqrt(k/12))
		}
	}
	else {
		print(method)
		stop('method is not recognized')
	}
	# return
	return(cmbd.pvalue)
	# gmeta.cmbd <- list(indiv.pvalues=gmi, method=method, cmbd.pvalue=c0mbd.pvalue)
	# return(gmeta.cmbd)
}
## [p-value combination tools]
### compute pdf of sum of n uniform(0,1), used in pvalue-based method
### input: q - quantile; n - number of uniform(0,1) r.v.s to sum.
### designed for n < 30, if n > 30 then use normal approximation.
psumunif <- function(q, n) {
	psunif <- 0
	for (i in 0:n) {
		psunif = psunif + (-1)^i * choose(n,i) * max(0, q - i)^n
	}
	psunif = psunif / factorial(n)
	psunif  
}
## [p-value combination tools - done]
## meta-analysis - combine p-value vector [done]



## meta-analysis - model based meta-analysis
# *****************************************************************************
#    gmeta.m() - model based meta-analysis methods
# *****************************************************************************
## main[gmeta.m]

## meta-analysis - model based meta-analysis
## main[gmeta.m]
gmeta.m <- function(gmi, gmi.type, method, linkfunc, weight, gmo.xgrid, ci.level, tau2, verbose) {
	gmeta.data <- gmeta.cdpvt.dproc(gmi, gmi.type, method, gmo.xgrid, tau2)
	gmeta.cmbd <- gmeta.cdpvt.combine(gmeta.data, method, linkfunc, weight, gmo.xgrid, ci.level, verbose)
	#class(gmeta.cmbd) <- c('gmeta.m')
	return(gmeta.cmbd)
}

### data processing
gmeta.cdpvt.dproc <- function(gmi, gmi.type, method, gmo.xgrid, tau2) {
	if (gmi.type == 'cd') {
		gmeta.data <- gmeta.cdpvt.dproc.cd(gmi, method, gmo.xgrid, tau2)
	} else { # gmi.type == 'pivot'
		gmeta.data <- gmeta.cdpvt.dproc.pivot(gmi, method, gmo.xgrid, tau2)
	}
	return(gmeta.data)
}
#####
gmeta.cdpvt.dproc.cd <- function(gmi, method, gmo.xgrid, tau2) {
	# unified data structure
	# x.grids and cdf values
	gmi.x     <- list()
	gmi.cd    <- list()
	# maintain mean & stddev
	nn1       <- length(gmi)
	gmi.theta <- numeric(nn1)
	gmi.sigma <- numeric(nn1)
	# unifying data structure:
	# x.grids, cdf-values-on-xgrids, mean, stddev for each study 
	for ( i in seq(1:nn1) ) {
		gmi.x[[i]]   <- gmi[[i]][,1]
		gmi.cd[[i]]  <- gmi[[i]][,2]
		gmi.theta[i] <- gmeta.cd.mean(gmi.x[[i]], gmi.cd[[i]])
		gmi.sigma[i] <- gmeta.cd.stddev(gmi.x[[i]], gmi.cd[[i]])
	}
	
	# calculate tau2
	if ( is.element(method, c('fixed-mle',
							  'fixed-robust1',
							  'fixed-robust2',
							  'fixed-robust2(sqrt12)')) ) {
		gmi.tau2 <- 0
	} else if ( is.element(method, c('random-mm',
									 'random-reml',
									#'random-tau2',
									 'random-robust1',
									 'random-robust2',
									 'random-robust2(sqrt12)')) ) {
		#gmi.tau2 <- Gtau2(gmi.theta, gmi.sigma, method)
		if ( is.numeric(tau2) ) {
			gmi.tau2 <- tau2
		} else if ( is.character(tau2) && is.element(tau2, c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB')) ) {
			gmi.tau2 <- Gtau2m(gmi.theta, gmi.sigma, tau2method=tau2)
		} else {
			gmi.tau2 <- Gtau2(gmi.theta, gmi.sigma, method)
		}
	} else if ( method == 'random-tau2' ) {
		if ( is.numeric(tau2) ) {
			gmi.tau2 <- tau2
		} else if ( is.character(tau2) && is.element(tau2, c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB')) ) {
			gmi.tau2 <- Gtau2m(gmi.theta, gmi.sigma, tau2method=tau2)
		} else {
			stop("tau2 not recongize: must specify a numeric tau2 or a tau2method in c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB').")
		}
	} else {
		stop('method not recongize.')
	}
	
	# unified x.grids
	# xl <- min(unlist(gmi.x))
	# xu <- max(unlist(gmi.x))
	# gmi.unix <- seq(from=xl, to=xu, length.out=n)
	gmi.unix <- gmo.xgrid # update in v2.0
	# calculate cdf-value on unified x.grids
	gmi.unix.cd <- NULL
	for ( i in seq(1,nn1) ) {
		gmi.unix.cd <- rbind(gmi.unix.cd, 
			approx(x=gmi.x[[i]], y=gmi.cd[[i]], xout=gmi.unix, rule=2)$y)
	}
	
	# return
	gmeta.data <- list( # x.grids
						gmi.x     = gmi.unix, # [unified x.grids - gmo.xgrids]
						# individual CDs
						gmi.cd    = gmi.unix.cd, 
						gmi.theta = gmi.theta, 
						gmi.sigma = gmi.sigma, 
						gmi.tau2  = gmi.tau2   )
	# return
	return(gmeta.data)
}
#####
gmeta.cdpvt.dproc.pivot <- function(gmi, method, gmo.xgrid, tau2) {
	# unified data structure
	# x.grids and cdf values
	#gmi.x     <- list()
	#gmi.cd    <- list()
	# maintain mean & stddev
	nn1       <- dim(gmi)[1]
	gmi.theta <- gmi[,1]
	gmi.sigma <- gmi[,2]
	
	# calculate tau2
	if ( is.element(method, c('fixed-mle',
							  'fixed-robust1',
							  'fixed-robust2',
							  'fixed-robust2(sqrt12)')) ) {
		gmi.tau2 <- 0
	} else if ( is.element(method, c('random-mm',
									 'random-reml',
									#'random-tau2',
									 'random-robust1',
									 'random-robust2',
									 'random-robust2(sqrt12)')) ) {
		#gmi.tau2 <- Gtau2(gmi.theta, gmi.sigma, method)
		if ( is.numeric(tau2) ) {
			gmi.tau2 <- tau2
		} else if ( is.character(tau2) && is.element(tau2, c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB')) ) {
			gmi.tau2 <- Gtau2m(gmi.theta, gmi.sigma, tau2method=tau2)
		} else {
			gmi.tau2 <- Gtau2(gmi.theta, gmi.sigma, method)
		}
	} else if ( method == 'random-tau2' ) {
		if ( is.null(tau2) ) {
			stop("tau2 must not NULL if method='random-tau2'")
		}
		if ( is.numeric(tau2) ) {
			gmi.tau2 <- tau2
		} else if ( is.character(tau2) && is.element(tau2, c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB')) ) {
			gmi.tau2 <- Gtau2m(gmi.theta, gmi.sigma, tau2method=tau2)
		} else {
			stop("tau2 not recongize: must specify a numeric tau2 or a tau2method in c('DL', 'HS', 'SJ', 'HE', 'ML', 'REML', 'EB').")
		}
	} else {
		stop('method not recongize.')
	}
	
	# unified x.grids
	# xl <- min(unlist(gmi.x))
	# xu <- max(unlist(gmi.x))
	# gmi.unix <- seq(from=xl, to=xu, length.out=n)
	gmi.unix <- gmo.xgrid # update in v2.0
	# calculate cdf-value on unified x.grids
	gmi.unix.cd <- NULL
	if (is.element(method, c('fixed-mle',
							 'fixed-robust1',
							 'fixed-robust2',
							 'fixed-robust2(sqrt12)'))) {
		for (i in seq(1,nn1)) {
			gmi.unix.cd <- rbind(gmi.unix.cd,
				pnorm(gmi.unix, gmi.theta[i], gmi.sigma[i]))
		}
	} else if (is.element(method, c('random-mm',
									'random-reml',
									'random-tau2',
									'random-robust1',
									'random-robust2',
									'random-robust2(sqrt12)'))) {
		for (i in seq(1,nn1)) {
			gmi.unix.cd <- rbind(gmi.unix.cd,
				pnorm(gmi.unix, gmi.theta[i], sqrt(gmi.sigma[i]^2+gmi.tau2)))
		}
	} else {
		stop('method not recongize.')
	}
	
	# return
	gmeta.data <- list( # x.grids
						gmi.x     = gmi.unix, # [unified x.grids - gmo.xgrids]
						# individual CDs
						gmi.cd    = gmi.unix.cd, 
						gmi.theta = gmi.theta, 
						gmi.sigma = gmi.sigma, 
						gmi.tau2  = gmi.tau2   )
	# return
	return(gmeta.data)
}

### combine: gaussian & double exponential
gmeta.cdpvt.combine <- function(gmeta.data, method, linkfunc, weight, gmo.xgrid, ci.level, verbose) {
	# use linkfunc combine
	if ( linkfunc == 'inverse-normal-cdf' ) {
		cmbdF <- gmeta.cdpvt.combine.gaussian(gmeta.data, method, weight, gmo.xgrid, verbose)
	} else if ( linkfunc == 'inverse-laplace-cdf' ) {
		cmbdF <- gmeta.cdpvt.combine.de(gmeta.data, method, weight, gmo.xgrid, verbose)
	} else {
		stop('linkfunc not recongize.')
	}
	
	# post-processing
	cmbdf   <- F2f(gmeta.data$gmi.x, cmbdF)
	cmbdmn  <- gmeta.cd.mean(gmeta.data$gmi.x, cmbdF)
	cmbdmdn <- gmeta.cd.median(gmeta.data$gmi.x, cmbdF)
	cmbdsdv <- gmeta.cd.stddev(gmeta.data$gmi.x, cmbdF)
	
	# post-processing on individual CDs
	# number of study
	K = dim(gmeta.data$gmi.cd)[1]
	# significance level
	alpha = 1 - ci.level
	# calculate individual and combined CIs
	indiv.ci <- NULL
	for (i in 1:K) {
		indiv.ci <- rbind(indiv.ci,
			gmeta.cd.mdncis(gmeta.data$gmi.x, gmeta.data$gmi.cd[i, ], alpha))
	}
	combined.ci <- gmeta.cd.mdncis(gmeta.data$gmi.x, cmbdF, alpha)[c(1,3)]
	
	# return
	gmeta.cmbd <- list( # x.grids
						x.grids            = gmeta.data$gmi.x,
						# individual CDs
						individual.cds     = gmeta.data$gmi.cd,
						individual.means   = gmeta.data$gmi.theta,
						individual.stddevs = gmeta.data$gmi.sigma,
						individual.medians = indiv.ci[,2], 
						individual.cis     = indiv.ci[,c(1,3)],
						# combined CD
						combined.cd        = cmbdF, 
						combined.density   = cmbdf, 
						combined.mean      = cmbdmn, 
						combined.sd        = cmbdsdv,
						combined.median    = cmbdmdn, 
						combined.ci        = combined.ci,
						# other information
						method             = method, 
						linkfunc           = linkfunc, 
						weight             = weight,
						tau2               = gmeta.data$gmi.tau2,
						ci.level           = ci.level,             
						verbose            = verbose              )
	# return
	return(gmeta.cmbd)
}
#####
gmeta.cdpvt.combine.gaussian <- function(gmeta.data, method, weight, gmo.xgrid, verbose) {
	# extract data structure
	gmi.x     <- gmeta.data$gmi.x
	gmi.cd    <- gmeta.data$gmi.cd
	gmi.tau2  <- gmeta.data$gmi.tau2
	gmi.theta <- gmeta.data$gmi.theta
	gmi.sigma <- gmeta.data$gmi.sigma
	
	# specify study-specific weight
	if ( is.null(weight) ) {
		if ( verbose ) {
			cat('\nlinkfunc is inverse-normal-cdf - default weight is inverse-variance.\n') #commentOut
		}
		weight <- 1 / sqrt(gmi.sigma^2 + gmi.tau2) # tau2 is zero if fixed-effect model, method suppress tau2 if random-effects model
	} else {
		if ( verbose ) {
			cat('\nuser-specified study-specific weights.\n') #commentOut
		} else {
			#cat('\nuser-specified study-specific weights.\n') #commentOut
		}
	}
	
	# meta via CD-framework
	if (method == 'fixed-mle') {
		#adptws <- 1 / gmi.sigma
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(weight*qnorm(gmi.cd), 2, sum) / sqrt(weight.ss))
	}
	#else if (method == 'fixed-bayesian') {
	#	#adptws <- 1 / gmi.sigma
	#	#weight <- weight * adptws [no adptws in v2.0]
	#	weight.ss <- sum(weight^2)
	#	cmbdF = pnorm(apply(weight*qnorm(gmi.cd), 2, sum) / sqrt(weight.ss))
	#}
	else if (method == 'fixed-robust1') {
		k <- length(weight)
		#adptws <- 1 / gmi.sigma
		#weight <- weight * adptws [no adptws in v2.0]
		e.kernel <- function(u00, h) { 
			1/h * dnorm(u00/h) 
		}
		data.iqr <- qnorm(0.75, gmi.theta, gmi.sigma) - qnorm(0.25, gmi.theta, gmi.sigma)
		cmbdm.all <- NULL
		cmbdF.all <- NULL
		for (i in 1:k) {
			# ith kernel weights
			knwsi <- e.kernel(gmi.theta-gmi.theta[i], sqrt(data.iqr))
			# adjusted weights [adjusted by kernel]
			wts.adj <- weight * knwsi
			weight.ss <- sum(wts.adj^2)
			# ith kernel combined Fs
			Fknl <- pnorm(apply(wts.adj*qnorm(gmi.cd), 2, sum)/sqrt(weight.ss))
			# ith kernel combined F's median
			mknl <- gmeta.cd.median(gmi.x, Fknl)
			# save for later analysis
			cmbdm.all <- c(cmbdm.all, mknl)
			cmbdF.all <- rbind(cmbdF.all, Fknl)
		}
		# using the median one as robust one
		srt.list <- sort.list(cmbdm.all)
		if (k %% 2) {
			# k is odd, median is (k+1)/2
			cmbdF <- cmbdF.all[srt.list[(k+1)/2], ]
		}
		else {
			# k is even, median is k/2 and k/2 + 1, combined with equal weight one
			cmbdF <- pnorm( (qnorm(cmbdF.all[srt.list[k/2], ]) + qnorm(cmbdF.all[srt.list[k/2+1], ]))/sqrt(2) )
		}
	}
	else if (method == 'fixed-robust2') {
		cmbdFr0 <- gmeta.cdpvt.combine.gaussian(gmeta.data, 
			method='fixed-robust2(sqrt12)', weight, gmo.xgrid, verbose)
		# latex format
		# $\hat{theta^{(c)}) = H(c)^{-1}(1/2)$
		x.median <- gmeta.cd.median(gmi.x, cmbdFr0)
		# ith cd at x.median
		nn1 <- length(weight)
		mdn.cds <- numeric(nn1)
		for (i in seq(1,nn1)) {
			mdn.cds[i] <- approx(x=gmi.x, y=gmi.cd[i,], xout=x.median)$y
		}
		weight.ss <- sum(weight^2*(mdn.cds-0.5)^2)
		cmbdF = pnorm(apply(weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'fixed-robust2(sqrt12)') {
		#adptws <- 1 / gmi.sigma
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(sqrt(12)*weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-mm') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(weight*qnorm(gmi.cd), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-reml') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(weight*qnorm(gmi.cd), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-tau2') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(weight*qnorm(gmi.cd), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-robust1') {
		k <- length(weight)
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		e.kernel <- function(u00, h) { 
			1/h * dnorm(u00/h) 
		}
		data.iqr <- qnorm(0.75, gmi.theta, gmi.sigma) - qnorm(0.25, gmi.theta, gmi.sigma)
		cmbdm.all <- NULL
		cmbdF.all <- NULL
		for (i in 1:k) {
			# ith kernel weights
			knwsi <- e.kernel(gmi.theta-gmi.theta[i], sqrt(data.iqr))
			# adjusted weights
			wts.adj <- weight * knwsi
			weight.ss <- sum(wts.adj^2)
			# ith kernel combined Fs
			Fknl <- pnorm(apply(wts.adj*qnorm(gmi.cd), 2, sum)/sqrt(weight.ss))
			# ith kernel combined F's median
			mknl <- gmeta.cd.median(gmi.x, Fknl)
			# save for later analysis
			cmbdm.all <- c(cmbdm.all, mknl)
			cmbdF.all <- rbind(cmbdF.all, Fknl)
		}
		# using the median one as robust one
		srt.list <- sort.list(cmbdm.all)
		if (k %% 2) {
			# k is odd, median is (k+1)/2
			cmbdF <- cmbdF.all[srt.list[(k+1)/2], ]
		}
		else {
			# k is even, median is k/2 and k/2 + 1, combined with equal weight one
			cmbdF <- pnorm( (qnorm(cmbdF.all[srt.list[k/2], ]) + qnorm(cmbdF.all[srt.list[k/2+1], ]))/sqrt(2) )
		}
	}
	else if (method == 'random-robust2') {
		cmbdFr0 <- gmeta.cdpvt.combine.gaussian(gmeta.data, 
			method='random-robust2(sqrt12)', weight, gmo.xgrid, verbose)
		# latex format
		# $\hat{theta^{(c)}) = H(c)^{-1}(1/2)$
		x.median <- gmeta.cd.median(gmi.x, cmbdFr0)
		# ith cd at x.median
		nn1 <- length(weight)
		mdn.cds <- numeric(nn1)
		for (i in seq(1,nn1)) {
			mdn.cds[i] <- approx(x=gmi.x, y=gmi.cd[i,], xout=x.median)$y
		}
		weight.ss <- sum(weight^2*(mdn.cds-0.5)^2)
		cmbdF = pnorm(apply(weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-robust2(sqrt12)') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(sqrt(12)*weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else {
		print(method)
		stop('method is not recongized.')
	}
	# return
	return(cmbdF)
}
#####
gmeta.cdpvt.combine.de <- function(gmeta.data, method, weight, gmo.xgrid, verbose) {
	# extract data structure
	gmi.x     <- gmeta.data$gmi.x
	gmi.cd    <- gmeta.data$gmi.cd
	gmi.tau2  <- gmeta.data$gmi.tau2
	gmi.theta <- gmeta.data$gmi.theta
	gmi.sigma <- gmeta.data$gmi.sigma
	
	# specify study-specific weight
	if ( is.null(weight) ) {
		if ( is.element(method, c('fixed-mle',
								  'random-mm',
								  'random-reml',
								  'random-tau2')) ) {
			if ( verbose ) {
				cat('\nlinkfunc is inverse-laplace-cdf - default weight is all one (Bahadur Efficiency).\n') #commentOut
			}
			weight <- rep(1, length(gmi.theta)) # tau2 is zero if fixed-effect model, method suppress tau2 if random-effects model
		} else if ( is.element(method, c('fixed-robust1',
										 'fixed-robust2',
										 'fixed-robust2(sqrt12)',
										 'random-robust1',
										 'random-robust2',
										 'random-robust2(sqrt12)')) ) {
			if ( verbose ) {
				cat('\nif use robust method, linkfunc is forced to use inverse-normal-cdf and default weight is inverse-variance.\n') #commentOut
			}
			weight <- 1 / sqrt(gmi.sigma^2 + gmi.tau2) # tau2 is zero if fixed-effect model, method suppress tau2 if random-effects model
		} else {
			stop('method not recongize.')
		}
	} else {
		if ( ( verbose ) && !( weight == 1 ) ) {
			cat("\nwarning: user-specified study-specific weights are not recommend if take linkfunc='inverse-laplace-cdf', the default weights with all ones give Bahadur Efficiency.\n") #commentOut
		} else {
			#cat("\nwarning: use user-specified weight which is not recommend.\n") #commentOut
		}
	}
	
	# meta via CD-framework
	# F0: double.exponential
	# linkfunc: inverse-laplace-cdf
	glog <- function(x) {
		ifelse(x<=0, -5000, log(x))
	}
	qde <- function(x) {
		ifelse(x>0.5, -glog(2*(1-x)), glog(2*x))
	}
	rde <- function(k) {
		rexp(k)*(2*rbinom(k,1,0.5)-1)
	}
	# meta via CD-framework
	if (method == 'fixed-mle') {
		#adptws <- 1 / gmi.sigma
		#weight <- weight * adptws [no adptws in v2.0]
		#p.cmbd.de <- convolution.sim(rde, weight) 
		p.cmbd.de <- .convolution.de(rde, weight, verbose) # [update in v2.0: .convolution.de()]
		cmbd.q = apply(weight*qde(gmi.cd), 2, sum)
		cmbdF <- p.cmbd.de(cmbd.q)
	}
	#else if (method == 'fixed-bayesian') {
	#	#adptws <- 1 / gmi.sigma
	#	#weight <- weight * adptws [no adptws in v2.0]
	#	#p.cmbd.de <- convolution.sim(rde, weight)
	#	p.cmbd.de <- .convolution.de(rde, weight, verbose) # [update in v2.0: .convolution.de()]
	#	cmbd.q = apply(weight*qde(gmi.cd), 2, sum)
	#	cmbdF <- p.cmbd.de(cmbd.q)
	#}
	else if (method == 'fixed-robust1') {
		k <- length(weight)
		#adptws <- 1 / gmi.sigma
		#weight <- weight * adptws [no adptws in v2.0]
		e.kernel <- function(u00, h) { 
			1/h * dnorm(u00/h) 
		}
		data.iqr <- qnorm(0.75, gmi.theta, gmi.sigma) - qnorm(0.25, gmi.theta, gmi.sigma)
		cmbdm.all <- NULL
		cmbdF.all <- NULL
		for (i in 1:k) {
			# ith kernel weights
			knwsi <- e.kernel(gmi.theta-gmi.theta[i], sqrt(data.iqr))
			# adjusted weights
			wts.adj <- weight * knwsi
			# ith kernel combined Fs
			p.cmbd.de <- convolution.sim(rde, wts.adj)
			#p.cmbd.de <- .convolution.de(rde, wts.adj, verbose=FALSE) # [update in v2.0: .convolution.de() - no need exact convolution since weight not all one]
			cmbd.q = apply(wts.adj*qde(gmi.cd), 2, sum)
			Fknl <- p.cmbd.de(cmbd.q)
			# ith kernel combined F's median
			mknl <- gmeta.cd.median(gmi.x, Fknl)
			# save for later analysis
			cmbdm.all <- c(cmbdm.all, mknl)
			cmbdF.all <- rbind(cmbdF.all, Fknl)
		}
		# using the median one as robust one
		srt.list <- sort.list(cmbdm.all)
		if (k %% 2) {
			# k is odd, median is (k+1)/2
			cmbdF <- cmbdF.all[srt.list[(k+1)/2], ]
		}
		else {
			# k is even, median is k/2 and k/2 + 1, combined with equal weight one
			cmbdF <- pnorm( (qnorm(cmbdF.all[srt.list[k/2], ]) + qnorm(cmbdF.all[srt.list[k/2+1], ]))/sqrt(2) )
		}
	}
	else if (method == 'fixed-robust2') {
		cmbdFr0 <- gmeta.cdpvt.combine.de(gmeta.data, 
			method='fixed-robust2(sqrt12)', weight, gmo.xgrid, verbose)
		# latex format
		# $\hat{theta^{(c)}) = H(c)^{-1}(1/2)$
		x.median <- gmeta.cd.median(gmi.x, cmbdFr0)
		# ith cd at x.median
		nn1 <- length(weight)
		mdn.cds <- numeric(nn1)
		for (i in seq(1,nn1)) {
			mdn.cds[i] <- approx(x=gmi.x, y=gmi.cd[i,], xout=x.median)$y
		}
		weight.ss <- sum(weight^2*(mdn.cds-0.5)^2)
		cmbdF = pnorm(apply(weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'fixed-robust2(sqrt12)') {
		#adptws <- 1 / gmi.sigma 
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(sqrt(12)*weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-mm') {
		#adptws <- 1 / sqrt((gmi.sigma^2 + gmi.tau2))
		#weight <- weight * adptws [no adptws in v2.0]
		#p.cmbd.de <- convolution.sim(rde, weight)
		p.cmbd.de <- .convolution.de(rde, weight, verbose) # [update in v2.0: .convolution.de()]
		cmbd.q = apply(weight*qde(gmi.cd), 2, sum)
		cmbdF <- p.cmbd.de(cmbd.q)
	}
	else if (method == 'random-reml') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		#p.cmbd.de <- convolution.sim(rde, weight)
		p.cmbd.de <- .convolution.de(rde, weight, verbose) # [update in v2.0: .convolution.de()]
		cmbd.q = apply(weight*qde(gmi.cd), 2, sum)
		cmbdF <- p.cmbd.de(cmbd.q)
	}
	else if (method == 'random-tau2') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		#p.cmbd.de <- convolution.sim(rde, weight)
		p.cmbd.de <- .convolution.de(rde, weight, verbose) # [update in v2.0: .convolution.de()]
		cmbd.q = apply(weight*qde(gmi.cd), 2, sum)
		cmbdF <- p.cmbd.de(cmbd.q)
	}
	else if (method == 'random-robust1') {
		k <- length(weight)
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		e.kernel <- function(u00, h) { 
			1/h * dnorm(u00/h) 
		}
		data.iqr <- qnorm(0.75, gmi.theta, gmi.sigma) - qnorm(0.25, gmi.theta, gmi.sigma)
		cmbdm.all <- cmbdF.all <- NULL
		for (i in 1:k) {
			# ith kernel weights
			knwsi <- e.kernel(gmi.theta-gmi.theta[i], sqrt(data.iqr))
			# adjusted weights
			wts.adj <- weight * knwsi
			# ith kernel combined Fs
			p.cmbd.de <- convolution.sim(rde, wts.adj)
			#p.cmbd.de <- .convolution.de(rde, wts.adj, verbose=FALSE) # [update in v2.0: .convolution.de() - no need exact convolution since weight not all one]
			cmbd.q = apply(wts.adj*qde(gmi.cd), 2, sum)
			Fknl <- p.cmbd.de(cmbd.q)
			# ith kernel combined F's median
			mknl <- gmeta.cd.median(gmi.x, Fknl)
			# save for later analysis
			cmbdm.all <- c(cmbdm.all, mknl)
			cmbdF.all <- rbind(cmbdF.all, Fknl)
		}
		# using the median one as robust one
		srt.list <- sort.list(cmbdm.all)
		if (k %% 2) {
			# k is odd, median is (k+1)/2
			cmbdF <- cmbdF.all[srt.list[(k+1)/2], ]
		}
		else {
			# k is even, median is k/2 and k/2 + 1, combined with equal weight one
			cmbdF <- pnorm( (qnorm(cmbdF.all[srt.list[k/2], ]) + qnorm(cmbdF.all[srt.list[k/2+1], ]))/sqrt(2) )
		}
	}
	else if (method == 'random-robust2') {
		cmbdFr0 <- gmeta.cdpvt.combine.de(gmeta.data, 
			method='random-robust2(sqrt12)', weight, gmo.xgrid, verbose)
		# latex format
		# $\hat{theta^{(c)}) = H(c)^{-1}(1/2)$
		x.median <- gmeta.cd.median(gmi.x, cmbdFr0)
		# ith cd at x.median
		nn1 <- length(weight)
		mdn.cds <- numeric(nn1)
		for (i in seq(1,nn1)) {
			mdn.cds[i] <- approx(x=gmi.x, y=gmi.cd[i,], xout=x.median)$y
		}
		weight.ss <- sum(weight^2*(mdn.cds-0.5)^2)
		cmbdF = pnorm(apply(weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else if (method == 'random-robust2(sqrt12)') {
		#adptws <- 1 / sqrt(gmi.sigma^2 + gmi.tau2)
		#weight <- weight * adptws [no adptws in v2.0]
		weight.ss <- sum(weight^2)
		cmbdF = pnorm(apply(sqrt(12)*weight*(gmi.cd-0.5), 2, sum)/sqrt(weight.ss))
	}
	else {
		print(method)
		stop('method is not recongized.')
	}
	# return
	return(cmbdF)
}
## meta-analysis - model based meta-analysis[done]

## meta-analysis - model based meta-analysis[done]



## meta-analysis - combine evidence from 2x2 tables
# *****************************************************************************
#    gmeta.e() - combine evidence from 2x2 tables (exact methods)
# *****************************************************************************
## main[gmeta.e]
gmeta.e <- function(gmi, method, weight, gmo.xgrid, ci.level, mc.iteration, eta, verbose, report.error) {
	data_matrix = as.matrix(cbind(gmi[,1], gmi[,3], gmi[,2], gmi[,4], gmi[,1]+gmi[,3]))
	colnames(data_matrix) = c("x","y","N","M","x+y")
	if (method == 'exact1') {
		gmeta.data <- gmeta.exact.indiv(data_matrix, gmo.xgrid, ci.level)
		gmeta.cmbd <- gmeta.exact.combine(gmeta.data, weight, gmo.xgrid, ci.level, mc.iteration, report.error)
	} else if (method == 'exact2') {
		gmeta.cmbd <- gmeta.exact.LT(data_matrix, weight, gmo.xgrid, ci.level, mc.iteration, eta, verbose, report.error)
	} else if (method == 'MH' || method == 'Mantel-Haenszel') {
		gmeta.cmbd <- gmeta.MH(data_matrix, weight, gmo.xgrid, ci.level)
	} else if (method == 'Peto') {
		gmeta.cmbd <- gmeta.peto(data_matrix, weight, gmo.xgrid, ci.level)
	} else {
		stop('gmi.type 2x2 only match methods MH, Mantel-Haenszel, Peto, exact1, exact2.')
	}
	#gmeta.cmbd$method <- method
	return(gmeta.cmbd)
}

### 2x2 with exact1(LLX)
##### data processing
gmeta.exact.indiv <- function(data_matrix, gmo.xgrid, ci.level) {
	# number of study
	K = nrow(data_matrix)
	# input  - individual study mean and standard deviation
	gmeta.theta = log( (data_matrix[,1]*(data_matrix[,4]-data_matrix[,2])) / (data_matrix[,2]*(data_matrix[,3]-data_matrix[,1])) ) 
	gmeta.sigma = sqrt(1/data_matrix[,1] + 1/data_matrix[,2] + 1/(data_matrix[,3]-data_matrix[,1]) + 1/(data_matrix[,4]-data_matrix[,2]))
	index = is.finite(gmeta.theta) & is.finite(gmeta.sigma)
	# output - individual study CDs
	indiv.cds = matrix(NA, K, length(gmo.xgrid))
	indiv.cis = matrix(NA, K, 2)
	indiv.medians = rep(NA, K)
	indiv.means = rep(NA, K)
	indiv.stddevs = rep(NA, K)
	# construct individual CDs
	alpha <- 1 - ci.level
	pfunc <- function(theta) {
		return(1 - pFNCHypergeo.wrap(data_matrix[i,1],data_matrix[i,3],data_matrix[i,4],data_matrix[i,5], theta) + 
		           dFNCHypergeo(data_matrix[i,1],data_matrix[i,3],data_matrix[i,4],data_matrix[i,5], theta) / 2)
	}
	for (i in 1:K) {
		indiv.cds[i, ] <- sapply(exp(gmo.xgrid), pfunc)
	}
	indiv.cds <- ifelse(indiv.cds < 1-0.1^12, indiv.cds, 1-0.1^12) # otherwise qnorm() will resolve Inf.
	# individual study - CI median
	for (i in 1:K) {
		indiv.cis[i, ]   <- log(c(.quantileCD(pfunc, alpha/2), .quantileCD(pfunc, 1-alpha/2)))
		indiv.medians[i] <- log(.quantileCD(pfunc, 0.5))
	}
	# individual study - mean and standard deviation
	na.moment <- is.na(gmeta.theta*0)
	for (i in c(1:K)[!na.moment]) {
		indiv.means[i]   <- gmeta.cd.mean(gmo.xgrid, indiv.cds[i,])
		indiv.stddevs[i] <- gmeta.cd.stddev(gmo.xgrid, indiv.cds[i,])
	}
	indiv.means[na.moment] <- gmeta.theta[na.moment]
	indiv.stddevs[na.moment] <- gmeta.sigma[na.moment]
	# return
	gmeta.data <- list( # input
						data_matrix   = data_matrix,
						# individual CDs
						indiv.cds     = indiv.cds,
						indiv.cis     = indiv.cis,
						indiv.medians = indiv.medians,
						indiv.means   = indiv.means,
						indiv.stddevs = indiv.stddevs,
						# output gridding points
						x.grids       = gmo.xgrid     )
	# return
	return(gmeta.data)
}
##### combine: exact1 method
gmeta.exact.combine <- function(gmeta.data, weight, gmo.xgrid, ci.level, mc.iteration, report.error) {
	#gmeta.cmbd <- gmeta.data
	#x.grids <- gmeta.data$x.grids
	# data_matrix
	data_matrix = gmeta.data$data_matrix
	# number of study
	K = nrow(data_matrix)
	# combine individual CDs
	alpha = 1 - ci.level
	xstmts <- .Estimates(data_matrix[,1:4])
	# obtain weight
	if ( is.null(weight) ) {
		weight <- xstmts$weight.hat
	} else {
		warning('combine 2x2 with method="exact1", default weight=NULL is strongly recommended.')
	}
	indiv.cds <- gmeta.data$indiv.cds
	# combined cdf
	combined.cd <- weight %*% qnorm( indiv.cds )
	combined.cd <- combined.cd / sqrt( sum( weight^2 ) )
	combined.cd <- as.vector( pnorm(combined.cd) ) # combined CD with values on limited points
	# combined CD function	
	combinedCDF <- function(theta) {  # combined CD function, used when searching for quantiles
		ff = sapply(1:K, FUN = function(j) { midp.oddsratio(data_matrix[j,1], data_matrix[j,3], data_matrix[j,4], data_matrix[j,5], theta) })
		ff = ifelse(ff<1-0.1^12, ff, 1-0.1^12)
		ff = pnorm( weight %*% qnorm(ff) / sqrt(sum(weight^2)) )
		return(ff)
	}
	# inference derived from combined CD function
	mn.xstmt = log(.quantileCD(combinedCDF, 1/2))            # mean
	lower.ci = log(.quantileCD(combinedCDF, alpha/2))        # ci.lower
	upper.ci = log(.quantileCD(combinedCDF, 1 - alpha/2))    # ci.upper
	o.pvalue = 2 * min( combinedCDF(1), 1 - combinedCDF(1) ) # null hypothesis: odd-ratio==1
	# adjust individual CDs and combined CD function to reduce test.size.error
	if( is.na(xstmts$or.hat)==TRUE || xstmts$or.hat==0 || xstmts$or.hat==Inf ) {
		mn.xstmt.adj = xstmts$or.hat
		lower.ci.adj = NA
		upper.ci.adj = NA
		o.pvalue.adj = NA
		coverage.prbblty.error     = NA
		coverage.prbblty.error.adj = NA
	} else {
		indiv.cds.adj   <- t( sapply(1:K, FUN = function(j) { adjust.beta(gmeta.data$indiv.cds[j,], data_matrix[j,3], xstmts$pt.hat[j], data_matrix[j,4], xstmts$pc.hat[j]) }) )
		combined.cd.adj <- weight %*% qnorm( indiv.cds.adj )
		combined.cd.adj <- combined.cd.adj / sqrt( sum( weight^2 ) )
		combined.cd.adj <- as.vector( pnorm(combined.cd.adj) ) # combined CD with values on limited points
		# combined CD function - adjusted
		combinedCDF.adj <- function(theta) { # combined CD function, used when searching for quantiles
			ff = sapply(1:K, FUN = function(j) { midp.oddsratio(data_matrix[j,1], data_matrix[j,3], data_matrix[j,4], data_matrix[j,5], theta) })
			ff = sapply(1:K, FUN = function(j) { adjust.beta(ff[j], data_matrix[j,3], xstmts$pt.hat[j], data_matrix[j,4], xstmts$pc.hat[j]) })
			ff = ifelse(ff<1-0.1^12, ff, 1-0.1^12)
			ff = pnorm( weight %*% qnorm(ff) / sqrt(sum(weight^2)) )
			return(ff)
		}
		# inference derived from combined CD function - adjusted
		mn.xstmt.adj = log(.quantileCD(combinedCDF.adj,1/2))               # mean
		lower.ci.adj = log(.quantileCD(combinedCDF.adj,alpha/2))           # ci.lower
		upper.ci.adj = log(.quantileCD(combinedCDF.adj,1-alpha/2))         # ci.upper
		o.pvalue.adj = 2 * min(combinedCDF.adj(1), 1 - combinedCDF.adj(1)) # null hypothesis: odd-ratio==1
		# calculate coverage probability error
		coverage.prbblty.error     <- 'Not Requested'
		coverage.prbblty.error.adj <- 'Not Requested'
		if ( report.error ) {	
			ee = test.size.error(alpha/2, data_matrix[,3], xstmts$pt.hat, data_matrix[,4], xstmts$pc.hat, weight, result.1minusS=TRUE, mc.iteration)
			coverage.prbblty.error     = ee$Test.Size.Error.1minusS - ee$Test.Size.Error
			coverage.prbblty.error.adj = ee$Test.Size.Error.Adjusted.1minusS - ee$Test.Size.Error.Adjusted
		}
	}
	# return
	gmeta.cmbd <- list( # input
						data_matrix = gmeta.data$data_matrix[,1:4],
						# individual CDs
						individual.cds     = gmeta.data$indiv.cds,
						individual.cis     = gmeta.data$indiv.cis,
						individual.medians = gmeta.data$indiv.medians,
						individual.means   = gmeta.data$indiv.means,
						individual.stddevs = gmeta.data$indiv.stddevs,
						# combined CD function
						combined.cd        = combined.cd,
						combined.density   = F2f(gmo.xgrid, combined.cd),
						combined.mean      = mn.xstmt,
						combined.median    = gmeta.cd.median(gmo.xgrid, combined.cd),
						combined.sd        = gmeta.cd.stddev(gmo.xgrid, combined.cd),
						combined.ci        = c(lower.ci, upper.ci),
						# p-value for null hypothesis: odd-ratio==1
						pvalue = o.pvalue,
						# combined CD function - adjusted 
						combined.cd.adjusted       = combined.cd.adj,
						combined.density.adjusted  = F2f(gmo.xgrid, combined.cd.adj),
						combined.mean.adjusted     = mn.xstmt.adj,
						combined.median.adjusted   = gmeta.cd.median(gmo.xgrid, combined.cd.adj),
						combined.sd.adjusted       = gmeta.cd.stddev(gmo.xgrid, combined.cd.adj),
						combined.ci.adjusted       = c(lower.ci.adj, upper.ci.adj),
						# p-value for null hypothesis: odd-ratio==1 - adjusted
						pvalue.adj = o.pvalue.adj,
						# coverage probability error - check liu2012exact
						coverage.prbblty.error     = coverage.prbblty.error,
						coverage.prbblty.error.adj = coverage.prbblty.error.adj,
						# other information
						method       = 'exact1',
						linkfunc     = 'inverse-fisher-exact-test-function', #[? adjusted]
						weight       = weight,
						tau2         = NULL,
						ci.level     = ci.level,
						verbose      = report.error,
						mc.iteration = mc.iteration,
						report.error = report.error,
						# output gridding points
						x.grids      = gmo.xgrid )
	# return
	return(gmeta.cmbd)
}
##### overall test size error for exact1 method
test.size.error <- function(s, n.vec, p1.vec, m.vec, p0.vec, weight.vec, result.1minusS=FALSE, mc.iteration=1000000) {
	# number of study
	K = length(n.vec)
	psi = (p1.vec[1]/(1-p1.vec[1])) / (p0.vec[1]/(1-p0.vec[1]))
	# generate p.i(\psi), i=1,2,...K. 
	mid.p.sample     = matrix(NA, nrow=mc.iteration, ncol=K) # each column is for an i
	mid.p.sample.adj = matrix(NA, nrow=mc.iteration, ncol=K) # each column is for an i
	for(i in 1:K){											
	    x = rbinom(mc.iteration, n.vec[i], p1.vec[i])
		y = rbinom(mc.iteration,m.vec[i],p0.vec[i])
		mid.p.sample[,i]     = midp.oddsratio(x, n.vec[i], m.vec[i], x + y, or=psi)
		mid.p.sample.adj[,i] = adjust.beta(mid.p.sample[,i], n=n.vec[i], pt=p1.vec[i], m=m.vec[i], pc=p0.vec[i])
	}
	# uniform sample
	uniform.sample = matrix(runif(mc.iteration*K,0,1), nrow=mc.iteration, ncol=K)
	# quantile
	qnorm.mid.p.sample     = qnorm(mid.p.sample)
	qnorm.mid.p.sample.adj = qnorm(mid.p.sample.adj)
	qnorm.uniform.sample   = qnorm(uniform.sample)
	# 
	bb1 = sapply(1:K, function(i) { ww<-weight.vec/weight.vec[i]; ww[1:K<=i]<-0; ww })
	bb2 = sapply(1:K, function(i) { ww<-weight.vec/weight.vec[i]; ww[1:K>=i]<-0; ww })
	#
	xoffset = sqrt(sum(weight.vec^2))/weight.vec*qnorm(s)
	#
	main     = qnorm.mid.p.sample     %*% bb1 + qnorm.uniform.sample %*% bb2
	main.adj = qnorm.mid.p.sample.adj %*% bb1 + qnorm.uniform.sample %*% bb2
	#
	aa     = t(xoffset-t(main))
	aa.adj = t(xoffset-t(main.adj))
	#
	pnorm.aa     = pnorm(aa)
	pnorm.aa.adj = pnorm(aa.adj)		
	#
	dd     = sapply(1:K, FUN = function(i) { findInterval(pnorm.aa[,i], sort(mid.p.sample[,i])) / mc.iteration - pnorm.aa[,i] })
	dd.adj = sapply(1:K, FUN = function(i) { findInterval(pnorm.aa.adj[,i], sort(mid.p.sample.adj[,i])) / mc.iteration - pnorm.aa.adj[,i] })
    # calculate test.size.error
	test.size.error     = sum(colMeans(dd))
	test.size.error.adj = sum(colMeans(dd.adj))
	# calculate test.size.error.1minusS
	test.size.error.1minusS     = 'Not requested'
	test.size.error.adj.1minusS = 'Not requested'
	if( result.1minusS ) {
		xoffset = sqrt(sum(weight.vec^2))/weight.vec*qnorm(1-s)  # use 1 - s instead of s
		#
		aa      = t(xoffset-t(main))
		aa.adj  = t(xoffset-t(main.adj))
		#
		pnorm.aa     = pnorm(aa)
		pnorm.aa.adj = pnorm(aa.adj)		
		#
		dd     = sapply(1:K, FUN = function(i) { findInterval(pnorm.aa[,i], sort(mid.p.sample[,i])) / mc.iteration - pnorm.aa[,i] })
		dd.adj = sapply(1:K, FUN = function(i) { findInterval(pnorm.aa.adj[,i], sort(mid.p.sample.adj[,i])) / mc.iteration - pnorm.aa.adj[,i] })
		# calculate test.size.error.1minusS
		test.size.error.1minusS     = sum(colMeans(dd))
		test.size.error.adj.1minusS = sum(colMeans(dd.adj))
	}
	# return
	return( list( Test.Size.Error                  = test.size.error, 
				  Test.Size.Error.Adjusted         = test.size.error.adj,
	              Test.Size.Error.1minusS          = test.size.error.1minusS,
				  Test.Size.Error.Adjusted.1minusS = test.size.error.adj.1minusS ) )
}

### 2x2 with exact2(LT's method) - risk difference[delta = p1 - p2]
gmeta.exact.LT <- function(data_matrix, weight, gmo.xgrid, ci.level, mc.iteration, eta, verbose, report.error) {
	# n=length(gmo.xgrid) - number of gridding delta [output]
	# mc.iteration - number of simulation for null distribution
	
	#if (!is.null(gmo.xgrid)) { warning('input gmo.xgrid will not be utilized, output range and griddings will be automatically set.') }#
	
	# expand data_matrix and set gmo.xgrid
	#data.mi = data_matrix[,1:4]
	data.mi = data_matrix[,c(2,1,4,3)] # [update in v2.0, target p1-p2, following program is on p2-p1, so change case-ctrl order]
	
	# number of trial
	nstudy = dim(data.mi)[1]
	
	# delta - individual studies
	n1 = data.mi[,3]
	n2 = data.mi[,4]
    p1 = data.mi[,1]/data.mi[,3]
	p2 = data.mi[,2]/data.mi[,4]
    deltap = p2 - p1
	
	# var weight - non-zero-event studies
	id = (1:nstudy)[p1*p2==0]
    n1[id] = data.mi[id,3]+1
	n2[id] = data.mi[id,4]+1
	p1[id] = (data.mi[id,1]+0.5)/(data.mi[id,3]+1)
	p2[id] = (data.mi[id,2]+0.5)/(data.mi[id,4]+1)
    varp = p1*(1-p1)/n1+p2*(1-p2)/n2
	wght = (n1*n2/(n1+n2))/sum(n1*n2/(n1+n2))
	
	# Mantel-Haenszel's method
	mu.MH = sum(deltap*wght)
	sd.MH = sqrt(sum(wght^2*varp))
    ci.MH = c(mu.MH-1.96*sd.MH, mu.MH+1.96*sd.MH)
    p.MH  = 1 - pchisq(mu.MH^2/sd.MH^2,1)
	
	# gridding detla
    d0 = max(abs(ci.MH))
	delta.grd = sort(c(0, seq(from=max(-1,-d0*15), to=min(1,d0*15), length=length(gmo.xgrid)-1)))
	
	# exact p-values [for observed data] given true delta[risk difference]
    diff.exact = function(x1, x2, n1, n2, delta.grd, n.grd=15, midp=TRUE) {
        # fit
		fit = binom.confint(x1, n1, 0.9995, methods='exact')
		# l,u
		l = fit$lower
		u = fit$upper
		# grd
        p1.grd = seq(l, u, length=n.grd)
		pnull1.tot = matrix(0, n1 + 1, n.grd)
        for( b in 1:n.grd ) { 
			p1 = p1.grd[b]
			pnull1.tot[,b] = dbinom(c(0:n1), n1, p1)
		}
		# df & sd
        dfnull = matrix(0, n1 + 1, n2 + 1)
		sdnull = matrix(0, n1 + 1, n2 + 1)
        for(i in 0:n1) {
			p1 = (i + 0.5) / (n1 + 1)
			p2 = (c(0:n2) + 0.5) / (n2 + 1)
            dfnull[i+1, ] = c(0:n2) / n2 - i / n1
            sdnull[i+1, ] = sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
        }
		# pv1 & pv2
        pv1 = numeric(0)
		pv2 = numeric(0)
        for( theta in delta.grd ) {
			p1 = (x1 + 0.5) / (n1 + 1)
			p2 = (x2 + 0.5) / (n2 + 1)
            tt  = (x2/n2-x1/n1-theta) / sqrt(p1*(1-p1)/n1+p2*(1-p2)/n2)
			# null hypothesis detla==theta
			tnull   = (dfnull-theta) / sdnull
            pvalue1 = rep(0, n.grd)
			pvalue2 = rep(0, n.grd)
			error = 1e-6
            for ( b in 1:n.grd ) { 
				p1 = p1.grd[b]
				p2 = p1 + theta
                if ( p2 >= 0 && p2 <= 1) {
					pnull1 = pnull1.tot[,b]
					pnull2 = dbinom(c(0:n2), n2, p2)
                    n1.adj = n1+1-max(c(1,(1:(n1+1))[cumsum(sort(pnull1))<error]))
                    n2.adj = n2+1-max(c(1,(1:(n2+1))[cumsum(sort(pnull2))<error]))
                    id1 = order(pnull1)
					id2 = order(pnull2)
                    id1 = (id1[(n1+1):1])[1:n1.adj]
					id2 = (id2[(n2+1):1])[1:n2.adj]
                    pnull = pnull1[id1] %*% t(pnull2[id2])
                    if ( midp==TRUE ) {
						pvalue1[b] = sum(pnull[tnull[id1,id2]>tt])+sum(pnull[tnull[id1,id2]==tt])*0.5
                        pvalue2[b] = sum(pnull[tnull[id1,id2]<tt])+sum(pnull[tnull[id1,id2]==tt])*0.5
                    } else {
                        pvalue1[b] = sum(pnull[tnull[id1,id2]>tt])+sum(pnull[tnull[id1,id2]==tt])
                        pvalue2[b] = sum(pnull[tnull[id1,id2]<tt])+sum(pnull[tnull[id1,id2]==tt])
                    }
                }
            }
            pv1=c(pv1, max(pvalue1)+(1-0.9995))
			pv2=c(pv2, max(pvalue2)+(1-0.9995))
        }
        return(list(pv1=pv1, pv2=pv2))
    }
	
	# pv1.pool&pv2.pool are the CDs of individual studies.
    pv1.pool = numeric(0)
	pv2.pool = numeric(0)
    for ( kk in 1:nstudy ) {
		# verbose
		if ( verbose ) {
			cat('2x2 exact2 processing trial:', kk, '\n') # use report.error as a surrogate as verbose [use parameter verbose - update v2.0]
		} # update status - due to slow speed of diff.exact().
		# resolve
		x1 = data.mi[kk,1]
		x2 = data.mi[kk,2]
		n1 = data.mi[kk,3]
		n2 = data.mi[kk,4]
		# fit
        fit = diff.exact(x1,x2,n1,n2,delta.grd,n.grd=15, midp=TRUE)
        pv1.pool = rbind(pv1.pool, fit$pv1)
		pv2.pool = rbind(pv2.pool, fit$pv2)
    }
	n = length(gmo.xgrid)
    for ( i in 1:nstudy ) {
		for ( j in 1:length(gmo.xgrid)) {
			pv1.pool[i,(n-j+1)] = max(pv1.pool[i,1:(n-j+1)])
			pv2.pool[i,j]       = max(pv2.pool[i,j:n])
        }
    }
	# pv1.pool&pv2.pool are the CDs of individual studies.
	
	# weight for individual study
	if ( is.null(weight) ) {
		weight = data.mi[,3] + data.mi[,4]
	} else {
		if ( verbose ) { # use report.error as a surrogate as verbose [use parameter verbose - update v2.0]
			cat('\nwarning: use user specified weights, instead of the default weight determined by study size.\n') 
		}
	}
	
	# combine individual CDs
	if ( is.numeric(eta) ) {
		# eta in [0,1]
		if ( ! (min(eta)>=0 && max(eta)<=1) ) {
			eta = seq(0.05, 0.95, length=20)
		}
		# F0^{-1}()
		F0inv <- function(ui) {
			sum( ((ui > (1-eta)) - eta) / ( eta*(1-eta) ) )
		}
	} else { # eta = 'Inf'
		# for the case K\to\inf:
		# pv1.pool&pv2.pool - adjust to make sense ln(ui/(1-ui))
		pv1.pool[pv1.pool<=0] = 1e-6
		pv1.pool[pv1.pool>=1] = 1 - 1e-6
		pv2.pool[pv2.pool<=0] = 1e-6
		pv2.pool[pv2.pool>=1] = 1 - 1e-6
		# pv2.pool = 1+(1-0.9995)*2 - pv2.pool
		# F0^{-1}()
		F0inv <- function(ui) {	log(ui/(1-ui)) }
	}

	# combine recipe g_c()
	gc <- function(u) { sum( sapply(u, F0inv) * weight ) }
	# simulation for G_c()
	T = numeric(mc.iteration)
	for (b in 1:mc.iteration) {
		u = runif(nstudy); T[b] = gc(u)
	}
	Gc = ecdf(T)
	
	# combined CDs
	Hc1 = Gc( apply(pv1.pool, 2, gc) ) # HC1 non-decreasing
	Hc2 = Gc( apply(pv2.pool, 2, gc) ) # HC2 non-increasing
	combined.ci = c(max(delta.grd[Hc1<0.025]), min(delta.grd[Hc2<0.025]))
	
	# post processing
	integrated2CDs <- function(cd1, cd2, delta.grd) {
		mdnpt = delta.grd[ round(median(which( abs(cd1-cd2) == min(abs(cd1-cd2)) ))) ]
		intgrtdcd = ifelse(delta.grd<mdnpt, cd1, 1+(1-0.9995)*2-cd2) # 1+(1-0.9995)*2-cd2
		for(j in 1:n) {
			intgrtdcd[j] = max(intgrtdcd[1:j])
		}
		intgrtdcd[intgrtdcd<=0] = 1e-6
		intgrtdcd[intgrtdcd>=1] = 1 - 1e-6;
		return(intgrtdcd)
	}
	
	# individual CDs
	indiv.cds <- NULL
	for (i in 1:nstudy) {
		indiv.cds = rbind(indiv.cds, integrated2CDs(pv1.pool[i,],pv2.pool[i,],delta.grd))
	}
	indiv.means   = numeric(0)
	indiv.stddevs = numeric(0)
	indiv.ci      = numeric(0)
	for (i in 1:nstudy) {
		indiv.means[i]   = gmeta.cd.mean(delta.grd, indiv.cds[i,])
		indiv.stddevs[i] = gmeta.cd.stddev(delta.grd, indiv.cds[i,])
		indiv.ci         = rbind(indiv.ci, gmeta.cd.mdncis(delta.grd, indiv.cds[i,], 1-ci.level))
	}
	
	# combined CD
	cmbdF   = integrated2CDs(Hc1, Hc2, delta.grd)
	cmbdf   = F2f(delta.grd, cmbdF)
	cmbdmn  = gmeta.cd.mean(delta.grd, cmbdF)
	cmbdmdn = gmeta.cd.median(delta.grd, cmbdF)
	cmbdsdv = gmeta.cd.stddev(delta.grd, cmbdF)
	# report.error
	if ( min(gmo.xgrid) < -1 || max(gmo.xgrid) > 1 ) {
		#gmo.xgrid = seq(from=-1, to=1, length=length(gmo.xgrid))
		gmo.xgrid = delta.grd
		if ( report.error ) {
			warning('\nuse exact2method combine 2x2 tables, parameter is risk-difference, only evaluate gmo.xgrid within [-1,1]\n')
		}
	}
	# return
	gmeta.cmbd <- list( # input
					    data_matrix = data.mi,
					    # individual CDs
					    individual.cds     = indiv.cds,
						individual.cis     = indiv.ci[,c(1,3)],
						individual.medians = indiv.ci[,2], 
					    individual.means   = indiv.means,
					    individual.stddevs = indiv.stddevs,
						# combined CD
						combined.cd        = cmbdF,
						combined.density   = cmbdf,
						combined.mean      = cmbdmn, 
						combined.median    = cmbdmdn,
						combined.sd        = cmbdsdv, 
						combined.ci        = combined.ci,		
						# other information
						method       = 'exact2',
						linkfunc     = 'log(u/(1-u))',
						weight       = weight,
						tau2         = NULL,
						ci.level     = ci.level,						
						mc.iteration = mc.iteration,
						eta          = eta,
						verbose      = verbose,
						report.error = report.error,
						# output gridding points
						x.grids      = delta.grd )
	# return
	return(gmeta.cmbd)
}

### Mantel-Haenszel's method
gmeta.MH <- function(data_matrix, weight, gmo.xgrid, ci.level) {
	# data matrix
	dm = data_matrix
	# number of study
	K = nrow(data_matrix)
	alpha = 1 - ci.level
	
	Ni = dm[,3]+dm[,4]
	
	R = dm[,1]*(dm[,4]-dm[,2])/Ni
	S = dm[,2]*(dm[,3]-dm[,1])/Ni
	P = (dm[,1]+dm[,4]-dm[,2])/Ni
	Q = (dm[,2]+dm[,3]-dm[,1])/Ni

	gmeta.theta = R/S
	gmeta.sigma = sqrt(1/dm[,1]+1/(dm[,3]-dm[,1])+1/dm[,2]+1/(dm[,4]-dm[,2])) * gmeta.theta
	
	# this is the same as var_U and var_US in Robins et.al(1986)'s paper when K=1.
	index = is.finite(gmeta.theta) & is.finite(gmeta.sigma)
	# individual CDs
	indiv.cds = matrix(NA,K,length(gmo.xgrid))
	indiv.cis = matrix(NA,K,2)
	for (ii in 1:K) {
		indiv.cds[ii,] = pnorm( ( gmo.xgrid - gmeta.theta[ii] ) / gmeta.sigma[ii] )
		indiv.cis[ii,] = c(gmeta.theta[ii]+qnorm(alpha/2)*gmeta.sigma[ii], gmeta.theta[ii]-qnorm(alpha/2)*gmeta.sigma[ii])
	}
	# combined CD
	cmbd.mean    = sum(R) / sum(S)
	cmbd.stddev  = sqrt(sum(P*R)/2/sum(R)^2+sum(P*S+Q*R)/2/sum(R)/sum(S)+sum(Q*S)/2/sum(S)^2) * cmbd.mean
	cmbd.cd      = pnorm( ( gmo.xgrid - cmbd.mean ) / cmbd.stddev )
	cmbd.density = dnorm( ( gmo.xgrid - cmbd.mean ) / cmbd.stddev ) / cmbd.stddev
	cmbd.ci      = c(cmbd.mean+qnorm(alpha/2)*cmbd.stddev, cmbd.mean-qnorm(alpha/2)*cmbd.stddev)
	# weight - MH default
	if (!is.null(weight)) { warning('user specific weight will be surrogated by MH default weight.') }
	weight = gmeta.sigma*S / (cmbd.stddev*sum(S))
	# return
	gmeta.cmbd <- list( # input
						data_matrix = data_matrix, 
						# individual CDs
						individual.cds     = indiv.cds,
						individual.cis     = indiv.cis,
						individual.medians = gmeta.theta,						
						individual.means   = gmeta.theta,
						individual.stddevs = gmeta.sigma,
						# combined CD
						combined.cd        = cmbd.cd,
						combined.density   = cmbd.density,
						combined.mean      = cmbd.mean,
						combined.median    = cmbd.mean,
						combined.sd        = cmbd.stddev,
						combined.ci        = cmbd.ci,
						# other information
						method       = 'Mantel-Haenszel',
						linkfunc     = 'inverse-normal-cdf',
						tau2         = NULL,
						weight       = weight,
						ci.level     = ci.level,
						verbose      = FALSE,
						# output gridding points
						x.grids      = gmo.xgrid  )
	# return
	return(gmeta.cmbd)
}

### Peto's method
gmeta.peto <- function(data_matrix, weight, gmo.xgrid, ci.level) {
	# data matrix
	dm = data_matrix
	# number of study
	K = nrow(data_matrix)
	alpha = 1 - ci.level
	
	Ni = dm[,3]+dm[,4]
	
	gmeta.theta = (dm[,1]-dm[,5]*dm[,3]/Ni)*Ni^2*(Ni-1) / (dm[,5]*dm[,3]*dm[,4]*(Ni-dm[,5]))
	gmeta.sigma = 1 / sqrt(dm[,5]*dm[,3]*dm[,4]*(Ni-dm[,5])/Ni^2/(Ni-1))
	
	index = is.finite(gmeta.theta) & is.finite(gmeta.sigma)
	# individual CDs
	indiv.cds = matrix(NA,K,length(gmo.xgrid))
	indiv.cis = matrix(NA,K,2)
	for (ii in 1:K) {
		indiv.cds[ii,] = pnorm((gmo.xgrid-gmeta.theta[ii])/gmeta.sigma[ii])
		indiv.cis[ii,] = c(gmeta.theta[ii]+qnorm(alpha/2)*gmeta.sigma[ii], gmeta.theta[ii]-qnorm(alpha/2)*gmeta.sigma[ii])
	}
	# combined CD
	cmbd.mean    = (sum(dm[,1])-sum(dm[,5]*dm[,3]/Ni))/sum(dm[,5]*dm[,3]*dm[,4]*(Ni-dm[,5])/Ni^2/(Ni-1))
	cmbd.stddev  = 1 / sqrt(sum(dm[,5]*dm[,3]*dm[,4]*(Ni-dm[,5])/Ni^2/(Ni-1)))
	cmbd.cd      = pnorm( ( gmo.xgrid - cmbd.mean ) / cmbd.stddev )
	cmbd.density = dnorm( ( gmo.xgrid - cmbd.mean ) / cmbd.stddev ) / cmbd.stddev
	cmbd.ci      = c(cmbd.mean+qnorm(alpha/2)*cmbd.stddev, cmbd.mean-qnorm(alpha/2)*cmbd.stddev)
	# weight - Peto's method default
	if (!is.null(weight)) { warning('user specific weight will be surrogated by Peto-s method default weight.') }
	weight = cmbd.stddev / gmeta.sigma
	# return 
	
	# return
	gmeta.cmbd <- list( # input
						data_matrix = data_matrix, 
						# individual CDs
						individual.cds     = indiv.cds,
						individual.cis     = indiv.cis,
						individual.medians = gmeta.theta,						
						individual.means   = gmeta.theta,
						individual.stddevs = gmeta.sigma,
						# combined CD
						combined.cd        = cmbd.cd,
						combined.density   = cmbd.density,
						combined.mean      = cmbd.mean,
						combined.median    = cmbd.mean,
						combined.sd        = cmbd.stddev,
						combined.ci        = cmbd.ci,
						# other information
						method       = 'Peto',
						linkfunc     = 'inverse-normal-cdf',
						weight       = weight,
						tau2         = NULL,
						ci.level     = ci.level,
						verbose      = FALSE,
						# output gridding points
						x.grids     = gmo.xgrid  )
	# return
	return(gmeta.cmbd)
}
## meta-analysis - combine evidence from 2x2 tables[done]

# unified meta-analysis[done]




# postprocessing

# *****************************************************************************
#    # postprocessing - print, summary, plot, etc..
# *****************************************************************************
# print, summary, plot, etc.

## print&summary [gmeta.p]
###print
print.gmeta.p <- function(x, ...) {
	# Title
	cat('\t\tP-value combination through CD-Framework\n')
	# Call
	cat('\nCall:\n')
	print(x$call) #[*$call]: type of language
	# Results
	cat('\nCombine Method:   ', x$method, '\n')
	cat('\nCombined p-value: ', x$combined.pvalue, '\n')
	# Details
	cat('\nIndividual p-values:\n')
	print(x$individual.pvalues) #[nicer visual display]
}
###summary
summary.gmeta.p <- function(object, ...) {
	# set
	object.sry <- object
	# set class
	class(object.sry) <- "summary.gmeta.p"
	# return
	return(object.sry)
}
###print of summary
print.summary.gmeta.p <- function(x, ...) {
	# Title
	cat('\t\tP-value combination through CD-Framework\n')
	# Call
	cat('\nCall:\n')
	print(x$call) #[*$call]: type of language
	# Results
	cat('\nCombine Method:   ', x$method, '\n')
	cat('\nCombined p-value: ', x$combined.pvalue, '\n')
	# Details
	cat('\nIndividual p-values:\n')
	print(x$individual.pvalues) #[nicer visual display]
}

## print&summary [gmeta.m]
###print
print.gmeta.m <- function(x, ...) {
	# Title
	cat('\t\tModel-Based Meta-Analysis through CD-Framework\n')
	# Call
	cat('\nCall:\n')
	print(x$call) #[*$call]: type of language
	# Results
	cat('\nSummary of Combined CD:\n')
	cmbd.cd.summary = data.frame( mean   = format(round(x$combined.mean,4),nsmall=4),
		                         median  = format(round(x$combined.median,4),nsmall=4),
		              standard.deviation = format(round(x$combined.sd,4),nsmall=4)    )
	row.names(cmbd.cd.summary) <- 'Combined CD'
	print(cmbd.cd.summary)
	# Details
	cat('\nCombined Confidence Distribution:\n')
	# construct combined CD
	cmbd.cd <- data.frame( x = x$x.grids, 
		density = x$combined.density, probability = x$combined.cd )
	# set lower/upper bound
	xl  <- min(x$x.grids)
	xu  <- max(x$x.grids)
	# count
	n.grids <- sum(cmbd.cd$x >= xl & cmbd.cd$x <= xu)
	# print head/tail 
	if (n.grids <= 10) {
		print(cmbd.cd) # simple print all x-cd-points
	} else {
		cmbd.cdhead <- cmbd.cd[1:5, ]
		cmbd.cdtail <- cmbd.cd[(n.grids-5):n.grids, ]
		names(cmbd.cdtail) <- c(' ', ' ', ' ') # avoid duplicate titles
		# print head/tail 5 x-cd-points
		# output
		print(cmbd.cdhead)
		cat('\t...\n\t...\n')
		print(cmbd.cdtail)
		#cat('\n          \n')
	}
}
###summary
summary.gmeta.m <- function(object, ...) {
	# set
	object.sry <- list()
	# set Call
	object.sry$call <- object$call
	# set combined CD
	object.sry$cmbd <- data.frame(  mean   = object$combined.mean, 
	                            median  = object$combined.median,
							    stddev  = object$combined.sd,
							   ci.lower = object$combined.ci[1],
							   ci.upper = object$combined.ci[2])
	row.names(object.sry$cmbd) <- 'Combined CD'
	# set individual CDs
	object.sry$idiv <- data.frame(  mean   = object$individual.means, 
	                            median  = object$individual.medians,
							    stddev  = object$individual.stddevs,
							   ci.lower = object$individual.cis[,1],
							   ci.upper = object$individual.cis[,2])
	# set individual study names
	row.names(object.sry$idiv) <- object$study.names
	# set ci.level
	object.sry$ci.level <- object$ci.level
	# set number of study
	object.sry$n.study <- dim(object.sry$idiv)[1]
	# set class
	class(object.sry) <- 'summary.gmeta.m'
	# return
	return(object.sry)
}
###print of summary
print.summary.gmeta.m <- function(x, ...) {
	# Title
	cat('\t\tModel-Based Meta-Analysis through CD-Framework\n')
	# Call
	cat('\nCall:\n')
	print(x$call) #[*$call]: type of language
	# Results
	cat('\nSummary of Combined CD:\n')
	print(x$cmbd)
	cat('\nConfidence level =', x$ci.level, '\n')
	# Details
	cat('\nSummary of Individual CDs:\n')
	print(x$idiv)
	cat('\nConfidence level =', x$ci.level, '\n')
}

## print&summary [gmeta.e]
###print
print.gmeta.e <- function(x, ...) {
	# Title
	cat('\t\tExact Meta-Analysis Approach through CD-Framework\n')
	# Call
	cat('\nCall:\n')
	print(x$call) #[*$call]: type of language
	# Results
	cat('\nSummary of Combined CD:\n')
	cmbd.cd.summary = data.frame( mean   = format(round(x$combined.mean,4),nsmall=4),
		                         median  = format(round(x$combined.median,4),nsmall=4),
		              standard.deviation = format(round(x$combined.sd,4),nsmall=4)    )
	row.names(cmbd.cd.summary) <- 'Combined CD'
	print(cmbd.cd.summary)
	# Details
	cat('\nCombined Confidence Distribution:\n')
	# construct combined CD
	cmbd.cd <- data.frame( x = x$x.grids, 
		density = x$combined.density, probability = x$combined.cd )
	# set lower/upper bound
	xl  <- min(x$x.grids)
	xu  <- max(x$x.grids)
	# count
	n.grids <- sum(cmbd.cd$x >= xl & cmbd.cd$x <= xu)
	# print head/tail 
	if (n.grids <= 10) {
		print(cmbd.cd) # simple print all x-cd-points
	} else {
		cmbd.cdhead <- cmbd.cd[1:5, ]
		cmbd.cdtail <- cmbd.cd[(n.grids-5):n.grids, ]
		names(cmbd.cdtail) <- c(' ', ' ', ' ') # avoid duplicate titles
		# print head/tail 5 x-cd-points
		# output
		print(cmbd.cdhead)
		cat('\t...\n\t...\n')
		print(cmbd.cdtail)
		#cat('\n          \n')
	}
}
###summary
summary.gmeta.e <- function(object, ...) {
	# set
	object.sry <- list()
	# set Call
	object.sry$call <- object$call
	# set combined CD
	object.sry$cmbd <- data.frame(  mean   = object$combined.mean, 
	                            median  = object$combined.median,
							    stddev  = object$combined.sd,
							   ci.lower = object$combined.ci[1],
							   ci.upper = object$combined.ci[2])
	row.names(object.sry$cmbd) <- 'Combined CD'
	# set individual CDs
	object.sry$idiv <- data.frame(  mean   = object$individual.means, 
	                            median  = object$individual.medians,
							    stddev  = object$individual.stddevs,
							   ci.lower = object$individual.cis[,1],
							   ci.upper = object$individual.cis[,2])							   
	# set individual study names
	row.names(object.sry$idiv) <- object$study.names
	# set ci.level
	object.sry$ci.level <- object$ci.level
	# set number of study
	object.sry$n.study <- dim(object.sry$idiv)[1]
	# set class
	class(object.sry) <- 'summary.gmeta.e'
	# return
	return(object.sry)
}
###print of summary
print.summary.gmeta.e <- function(x, ...) {
	# Title
	cat('\t\tExact Meta-Analysis Approach through CD-Framework\n')
	# Call
	cat('\nCall:\n')
	print(x$call) #[*$call]: type of language
	# Results
	cat('\nSummary of Combined CD:\n')
	print(x$cmbd)
	cat('\nConfidence level =', x$ci.level, '\n')
	# Details
	cat('\nSummary of Individual CDs:\n')
	print(x$idiv)
	cat('\nConfidence level =', x$ci.level, '\n')
}

## plot functions
plot.gmeta <- function(gmo, studies=NULL, 
	plot.option=c('confidence-density',
		'confidence-curve', 'cv',
		'confidence-distribution', 'cdf'), 
	type='l', xlab='x', ylab='density', xlim=NULL, ylim=NULL, ...) {
	UseMethod("plot", gmo);
}
	
###plot functions - model based meta-anlaysis
plot.gmeta.m <- function(x, studies=NULL, 
	plot.option=c('confidence-density',
				  'confidence-curve', 'cv',
				  'confidence-distribution', 'cdf'), 
	type='l', xlab='x', ylab='density', xlim=NULL, ylim=NULL, ...) {
	# match plot.option
	mfplot <- match.call()
	plot.option = match.arg(plot.option)
	
	# plot.individual.studies
	# if ( is.null(studies) ) {
		# # take all studies
		# gmi     <- x$input
		# nn1     <- ifelse(is.list(gmi), length(gmi), dim(gmi)[1])
		# studies <- c(1:nn1) 
	# }
	
	# plot via plot.option
	if ( plot.option == 'confidence-density' ) {
		gmeta.plot.cdd(x, studies, type, xlab, ylab, xlim, ylim, ...)
	} else if ( plot.option == 'cv' || plot.option == 'confidence-curve' ) {
		if ( ylab == 'density' ) {
			ylab = 'confidence curve'
		}
		gmeta.plot.cvs(x, studies, type, xlab, ylab, xlim, ylim, ...)
	} else if ( plot.option == 'cdf' || plot.option == 'confidence-distribution' ) {
		if ( ylab == 'density' ) {
			ylab = 'confidence distribution'
		}
		gmeta.plot.cdf(x, studies, type, xlab, ylab, xlim, ylim, ...)
	} else {
		stop('plot.option not recognize')
	}
}

###plot functions - combine evidence from 2x2 tables
plot.gmeta.e <- function(x, studies=NULL, 
	plot.option=c('confidence-density',
				  'confidence-curve', 'cv',
				  'confidence-distribution', 'cdf'), 
	type='l', xlab='x', ylab='confidence density', xlim=NULL, ylim=NULL, ...) {
	# match plot.option
	mfplot <- match.call()
	plot.option = match.arg(plot.option)
	
	# plot.individual.studies
	# if ( is.null(studies) ) {
		# # take all studies
		# gmi     <- x$input
		# nn1     <- ifelse(is.list(gmi), length(gmi), dim(gmi)[1])
		# studies <- c(1:nn1) 
	# }
	# take only non-zero/zero-event studies
	idx     <- studies[is.na(x$individual.mean)[studies]]
	studies <- studies[!is.na(x$individual.mean)[studies]]
	if ( length(idx) != 0 ) {
		for (i in idx) {
			cat('The confidence distribution of study', i, 'cannot be plot because it contains zero-zero events.','\n')
		}
	}
	# plot via plot.option
	if ( plot.option == 'confidence-density' ) {
		gmeta.plot.cdd(x, studies, type, xlab, ylab, xlim, ylim, ...)
	} else if ( plot.option == 'cv' || plot.option == 'confidence-curve' ) {
		if ( ylab == 'confidence density' ) {
			ylab = 'confidence curve'
		}
		gmeta.plot.cvs(x, studies, type, xlab, ylab, xlim, ylim, ...)
	} else if ( plot.option == 'cdf' || plot.option == 'confidence-distribution' ) {
		if ( ylab == 'confidence density' ) {
			ylab = 'confidence distribution'
		}
		gmeta.plot.cdf(x, studies, type, xlab, ylab, xlim, ylim, ...)
	} else {
		stop('plot.option not recognize')
	}
}

###plot functions - shared by model based and exact methods on 2x2 tables
#####
gmeta.plot.cdd <- function(x, studies, type, xlab, ylab, xlim, ylim, ...) {
	#e <- studies
	# number of layers
	nn <- length(studies) + 1 # individual studies + combined CD
	# extract study names
	enames <- c(x$study.names[studies], 'combined.density')
	# x.grids
	x.grids <- x$x.grids
	# individual CDs
	ecds <- x$individual.cds[studies, ]
	# individual densities
	edns <- NULL
	for (i in studies) {
		edns <- rbind(edns, F2f(x.grids, x$individual.cds[i, ]))
	}
	# individual and combined density
	edns <- rbind(edns, x$combined.density)
	# unified y-range on all layers for visual comparison
	if ( !is.null(studies) ) {
		if ( max(edns, na.rm=T) > 1 ) {
			edns <- edns / max(edns, na.rm=T)
		}
	}
	# individual and combined medians
	mdn <- c(x$individual.medians[studies], x$combined.median)
	# individual and combined confidence intervals
	emdncis <- rbind(x$individual.cis[studies,], x$combined.ci)
	# set lower/upper bound
	if ( !is.null(xlim) ) {
		xl <- min(xlim)
		xu <- max(xlim)
	} else {
		xl <- min(x.grids)
		xu <- max(x.grids)
	}
	# set study.names position
	xnames <- xl
	# plot
	if ( !is.null(studies) ) {
		plot(0, type='n', xlab=xlab, ylab=ylab, xlim=c(xl,xu), ylim=c(0,nn), yaxt='n')
	} else {
		plot(0, type='n', xlab=xlab, ylab=ylab, xlim=c(xl,xu), ylim=c(0,ceiling(max(edns))), yaxt='n')
	}
	for ( i in 1:nn ) {
		# supporting line
		abline(h=i-1, lty=2)
		# plot of density
		lines(x.grids, edns[i, ]+i-1, lty=1)
		# mark median, ci.lower, ci.upper
		points(mdn[i], i-1, cex=1, col='dark red')
		points(emdncis[i,1], i-1, pch='[', cex=1, col='dark red')
		points(emdncis[i,2], i-1, pch=']', cex=1, col='dark red')
		# legend study.names
		legend(xnames, i-0.5, legend=paste('#', enames[i], sep=''), bty='n')
	}
	# overall median
	abline(v=x$combined.median, lty=2)
}
#####
gmeta.plot.cvs <- function(gmo, studies, type, xlab, ylab, xlim, ylim, ...) {
	# number of layers
	nn <- length(studies) + 1 # individual studies + combined CD
	# extract study names
	enames <- c(gmo$study.names[studies], 'combined.cv')
	# x.grids
	x.grids <- gmo$x.grids
	
	# construct individual and combined CVs
	combined.cv    <- 1 - 2 * abs( gmo$combined.cd - 0.5 )
	individual.cvs <- 1 - 2 * abs( gmo$individual.cds - 0.5 )
	# set cvs
	ecvs <- rbind(individual.cvs[studies,], combined.cv)
	
	# individual and combined medians
	mdn <- c(gmo$individual.medians[studies], gmo$combined.median)
	# individual and combined confidence intervals
	emdncis <- rbind(gmo$individual.cis[studies,], gmo$combined.ci)
	# set lower/upper bound
	if ( !is.null(xlim) ) {
		xl <- min(xlim)
		xu <- max(xlim)
	} else {
		xl <- min(x.grids)
		xu <- max(x.grids)
	}
	# set study.names position
	xnames <- xl
	# plot
	if ( !is.null(studies) ) {
		plot(0, type='n', xlab=xlab, ylab=ylab, xlim=c(xl,xu), ylim=c(0,nn), yaxt='n')
	} else {
		plot(0, type='n', xlab=xlab, ylab=ylab, xlim=c(xl,xu), ylim=c(0,ceiling(max(ecvs))), yaxt='n')
	}
	for ( i in 1:nn ) {
		# supporting line
		abline(h=i-1, lty=2)
		# plot of density	
		lines(x.grids, ecvs[i, ]+i-1, lty=1)
		# mark median, ci.lower, ci.upper
		points(mdn[i], i-1, cex=1, col='dark red')
		points(emdncis[i,1], i-1, pch='[', cex=1, col='dark red')
		points(emdncis[i,2], i-1, pch=']', cex=1, col='dark red')
		# legend study.names
		legend(xnames, i-0.5, legend=paste('#', enames[i], sep=''), bty='n')
	}
	# overall median
	abline(v=gmo$combined.median, lty=2)
}
#####
gmeta.plot.cdf <- function(gmo, studies, type, xlab, ylab, xlim, ylim, ...) {
	#e <- studies
	# number of layers
	nn <- length(studies) + 1 # individual studies + combined CD
	# extract study names
	enames <- c(gmo$study.names[studies], 'combined.cdf')
	# x.grids
	x.grids <- gmo$x.grids
	# individual CDs
	ecds <- gmo$individual.cds[studies, ]
	# individual and combined CDs
	ecds <- rbind(ecds, gmo$combined.cd)
	# unified y-range on all layers for visual comparison
	if ( !is.null(studies) ) {
		if ( max(ecds, na.rm=T) > 1 ) {
			ecds <- ecds / max(ecds, na.rm=T)
		}
	}
	# individual and combined medians
	mdn <- c(gmo$individual.medians[studies], gmo$combined.median)
	# individual and combined confidence intervals
	emdncis <- rbind(gmo$individual.cis[studies,], gmo$combined.ci)
	# set lower/upper bound
	if ( !is.null(xlim) ) {
		xl <- min(xlim)
		xu <- max(xlim)
	} else {
		xl <- min(x.grids)
		xu <- max(x.grids)
	}
	# set study.names position
	xnames <- xl
	# plot
	if ( !is.null(studies) ) {
		plot(0, type='n', xlab=xlab, ylab=ylab, xlim=c(xl,xu), ylim=c(0,nn), yaxt='n')
	} else {
		plot(0, type='n', xlab=xlab, ylab=ylab, xlim=c(xl,xu), ylim=c(0,ceiling(max(ecds))), yaxt='n')
	}
	for ( i in 1:nn ) {
		# supporting line
		abline(h=i-1, lty=2)
		# plot of CDs
		lines(x.grids, ecds[i, ]+i-1, lty=1)
		# mark median, ci.lower, ci.upper
		points(mdn[i], i-1, cex=1, col='dark red')
		points(emdncis[i,1], i-1, pch='[', cex=1, col='dark red')
		points(emdncis[i,2], i-1, pch=']', cex=1, col='dark red')
		# legend study.names
		legend(xnames, i-0.5, legend=paste('#', enames[i], sep=''), bty='n')
	}
	# overall median
	abline(v=gmo$combined.median, lty=2)
}

# postprocessing [done]




# *****************************************************************************
#    other functions used in this package
# *****************************************************************************
# other functions used in this package

### using CD to make inference
### get mean, median, sd, ci, given cd, used everywhere
##### get mean from cd
gmeta.cd.mean <- function(x, cd) {
	x  =  x[!is.na(cd)]
	cd = cd[!is.na(cd)]	
	if ( length(unique(cd)) == 1 ) { 
		mn <- NA
	} else if (length(x) == length(cd)) {
		dn <- F2f(x,cd) # density
		mn <- sum(x*dn, na.rm=T) / sum(dn, na.rm=T) # mean
	} else {
		stop('length of x is not the same as length of cd.')
	}
	return(mn)
}
##### get median from cd
gmeta.cd.median <- function(x, cd) {
	x  =  x[!is.na(cd)]
	cd = cd[!is.na(cd)]	
	if ( length(unique(cd)) == 1 ) { 
		mdn <- NA
	} else {
		# same as approx in 
		# gmeta.cd.mdncis()
		x1 <- rev(x[cd<0.5])[1]
		x2 <- x[cd>0.5][1]
		y1 <- rev(cd[cd<0.5])[1]
		y2 <- cd[cd>0.5][1]
		# linear interpolation
		mdn <- x1+(x2-x1)/(y2-y1)*(0.5-y1)
	}
	return(mdn)
}
##### get stddev from cd
gmeta.cd.stddev <- function(x, cd) {
	x  =  x[!is.na(cd)]
	cd = cd[!is.na(cd)]	
	if ( length(unique(cd)) == 1 ) { 
		sdv <- NA
	} else {
		sdv <- diff(approx(x=cd, y=x, xout=c(0.25,0.75), ties='mean')$y) / (qnorm(0.75) - qnorm(0.25))
	}
	return(sdv)
}
##### get ci&median from cd
gmeta.cd.mdncis <- function(x, cd, alpha) {
	x  =  x[!is.na(cd)]
	cd = cd[!is.na(cd)]	
	if ( length(unique(cd)) == 1 ) { 
		mdncis <- c(NA, NA, NA)
	} else {
		mdncis <- approx(x=cd, y=x, xout=c(alpha/2, 0.5 ,1-alpha/2), ties='mean')$y
	}
	return(mdncis)
}

### computing density(pdf) from distribution(cdf), used everywhere
### mainly used for getting CD density from CD distribution, NA on two end to avoid sudden jump
F2f<-function (x, Fx) {
    fx = diff(Fx, lag = 2)/diff(x, lag = 2)
	fx = c(NA, fx, NA) # match the length of x-fx-cdf
	# n.xgrids
	n = length(fx)
	# liner extrapolation
	fx[1] = max(0, (fx[2]-fx[3])/(x[2]-x[3]) * (x[1]-x[2]) + fx[2])
	fx[n] = max(0, (fx[n-1]-fx[n-2])/(x[n-1]-x[n-2]) * (x[n]-x[n-1]) + fx[n-1])
	# return
	return(fx)
	#[update in v2.0: previously #return(c(NA, fx, NA)) # NA so that the edge of the plot of density is right.]
}



## [p-value combination tools]
## [update: p-value combination is implemented in C under this version]
## [p-value combination tools - done]



## [model based methods tools]

### computing tau2
### used in model based meta-analysis, random-effects models, estimate of tau2
##### Gtau2 - robust heterogeneity tau2 function
Gtau2 <- function(theta, sigma, method) {
	# calculate tau2 by method
	if ( method == 'random-mm' ) {
		tau2 <- GtauDL2(theta, sigma)	
	}
	else if ( method == 'random-reml' ) {
		tau2 <- GtauREML2(theta, sigma)
	}
	else if ( method == 'random-robust1' ) {
		tau2 <- GtauREML2(theta, sigma)	
	}
	else if ( method == 'random-robust2' ) {
		tau2 <- GtauREML2(theta, sigma)	
	}
	else if ( method == 'random-robust2(sqrt12)' ) {
		tau2 <- GtauREML2(theta, sigma)	
	}
	else {
		tau2 <- 0 # for fixed effects models
		#cat('\nmethod not recongize\n') #commentOut
	}
	# return
	return(tau2)
}
##### Gtau2m() methods for tau2
Gtau2m <- function(theta, sigma, tau2method) {
	# calculate tau2 by method-name
	if ( tau2method == 'DL' ) {
		tau2 <- GtauDL2(theta, sigma)
	}
	else if ( tau2method == 'HS' ) {
		tau2 <- GtauHS2(theta, sigma)
	}
	else if ( tau2method == 'SJ' ) {
		tau2 <- GtauSJ2(theta, sigma)
	}
	else if ( tau2method == 'HE' ) {
		tau2 <- GtauHE2(theta, sigma)
	}
	else if ( tau2method == 'ML' ) {
		tau2 <- GtauML2(theta, sigma)
	}
	else if ( tau2method == 'REML' ) {
		tau2 <- GtauREML2(theta, sigma)
	}
	else if ( tau2method == 'EB' ) {
		tau2 <- GtauEB2(theta, sigma)
	}
	else {
		stop('tau2method not recongize.')
	}
	# return
	return(tau2)
}
##### Get tauDL2
# ?update: need robust version
GtauDL2 <- function(theta, sigma) {
	nn1 = length(theta)
	wis <- 1 / sigma^2
	theta.mle <- sum(theta * wis) / sum(wis)
	q.w <- sum(wis * (theta - theta.mle)^2)
	tauDL2 <- max(0, (q.w - (nn1-1) ) / (sum(wis) - sum(wis^2)/sum(wis)) )
	return(tauDL2)
}
##### Get tauHS2
GtauHS2 <- function(theta, sigma) {
	nn1 = length(theta)
	wis <- 1 / sigma^2
	theta.mle <- sum(theta * wis) / sum(wis)
	q.w <- sum(wis * (theta - theta.mle)^2)
	tauHS2 <- max(0, (q.w - nn1) / sum(wis))
	return(tauHS2)
}
##### Get tauSJ2
GtauSJ2 <- function(theta, sigma) {
	nn1 = length(theta)
	tauSJ2 <- var(theta) * (nn1 - 1)/nn1
	wis <- 1 / (sigma^2 + tauSJ2)
	theta.mle <- sum(theta * wis) / sum(wis)
	q.w <- sum(wis * (theta - theta.mle)^2)
	tauSJ2 <- max(0, tauSJ2 * q.w / (nn1 - 1))
	return(tauSJ2)
}
##### Get tauHE2
GtauHE2 <- function(theta, sigma) {
	nn1 = length(theta)
	wis <- rep(1, nn1)
	theta.mle <- sum(theta * wis) / sum(wis)
	q.w <- sum(wis * (theta - theta.mle)^2)
	adj.w <- sum(sigma^2) * (nn1 - 1) / nn1
	tauHE2 <- max(0, (q.w - adj.w) / (nn1 - 1))
	return(tauHE2)
}
##### Get tauML2
GtauML2 <- function(theta, sigma) {
	nn1 = length(theta)
	# initialization
	tauML2 = max(0, GtauHE2(theta, sigma))
	# update
	nloop = 0
	absch = 1 # absolute change in tauML2 value
	while ( absch > 10^(-5) ) {
		nloop = nloop + 1
		if ( nloop > 10^5 ) {
			stop("tauML2 via MLE method does not converge.")
		} 
		else {
			tauML2O <- tauML2 # tauML2Old
			# update wML
			wML <- 1 / (sigma^2 + tauML2)
			WML <- diag(wML) - wML %o% wML / sum(wML)
			# update tauML2
			adj.tauML2 <- (t(theta) %*% WML %*% WML %*% theta - sum(wML)) / sum(wML^2)
			while ( tauML2 + adj.tauML2 < 0 ) {
				adj.tauML2 <- adj.tauML2 / 2
			}
			tauML2 <- tauML2 + adj.tauML2
			absch <- abs(tauML2 - tauML2O)
		}
	}
	# return
	return(max(0, tauML2))
}
##### Get tauREML2 
GtauREML2 <- function(theta, sigma) {
	nn1 = length(theta)
 	# adopt tauHE2 as initial tauREML2
	tauR2 <- GtauHE2(theta, sigma)
	# update
	nloop = 0
	absch = 1 # absolute change in tauR2 value
	while ( absch > 10^(-5) ) {
		nloop = nloop + 1
		if ( nloop > 10^5 ) {
			stop("tauREML2 via REML method does not converge.")
		}
		else {
			tauR2O <- tauR2 # tauR2Old
			# update thetaR, wR
			wR <- 1/(sigma^2 + tauR2O)
			thetaR <- sum(wR*theta) / sum(wR)
			# update tauR
			tauR2 <- sum(wR^2*(nn1/(nn1-1)*(theta-thetaR)^2 - sigma^2)) / sum(wR^2)
			absch <- abs(tauR2 - tauR2O)
		}
	}
	# return
	return(max(0, tauR2))
}
##### Get tauEB2
GtauEB2 <- function(theta, sigma) {
	nn1 = length(theta)
	# initialization
	tauEB2 <- max(0, GtauHE2(theta, sigma))
	# update
	nloop = 0
	absch = 1 # absolute change in tauEB2 value
	while ( absch > 10^(-5) ) {
		nloop = nloop + 1
		if ( nloop > 10^5 ) {
			stop("tauEB2 via EB method does not converge.")
		}
		else {
			tauEB2O <- tauEB2 # tauEB2Old
			# update thetaEB, wEB
			wEB <- 1 / (sigma^2 + tauEB2)
			thetaEB <- sum(theta * wEB) / sum(wEB)
			adj.w <- sum(wEB * (theta - thetaEB)^2)
			# update tauEB2
			adj.tauEB2 <- 1/sum(wEB)*(adj.w*nn1/(nn1-1)-nn1)
			while (tauEB2 + adj.tauEB2 < 0) {
				adj.tauEB2 <- adj.tauEB2 / 2
			}
			tauEB2 <- tauEB2 + adj.tauEB2
			absch <- abs(tauEB2 - tauEB2O)
		}
	}
	# return
	return(max(0, tauEB2))
}

### compute convolution distribution, used in double exp
convolution.sim <- function(rF0, weight, m=1e5) {
	k = length(weight)
	yy = numeric(m)
	for (i in c(1:m)) {
		yy[i] <- sum(rF0(k)*weight)
	}
	return(ecdf(yy))
}
### compute convolution distribution
### update in v2.0 - exact distribution functions
.convolution.de <- function(rF0, weight, verbose, m=1e5) {
	k = length(weight)
	if ( k <= 41 && all( weight == 1 ) ) {
		# return
		convolutionDEfunc <- .convolution.de.vfunc(k)
	} else {
		# verbose
		if ( verbose ) {
			if ( k > 41 ) {
				cat('\nlarge number of study, convolution distribution obtained by simulation.\n') #commentOut
			}
			if ( ! all( weight == 1 ) ) {
				cat('\nwarning: NOT all one weight, convolution distribution obtained by simulation.\n') #commentOut
			}
		}
		# return
		convolutionDEfunc <- convolution.sim(rF0, weight, m=m)
	}
	# return
	return( convolutionDEfunc )
}
# exact DE convolution disribution function k=1-41
.convolution.de.vfunc <- function(k) {
	# v1
	v1 <- function(x) {
		return( 1 )
	}
	# v2
	v2 <- function(x) {
		return( x/2 + 1 )
	}
	# v3
	v3 <- function(x) {
		return( x^2/8 + (5*x)/8 + 1 )
	}
	# v4
	v4 <- function(x) {
		return( x^3/48 + (3*x^2)/16 + (11*x)/16 + 1 )
	}
	# v5
	v5 <- function(x) {
		return( x^4/384 + (7*x^3)/192 + (29*x^2)/128 + (93*x)/128 + 1 )
	}
	# v6
	v6 <- function(x) {
		return( x^5/3840 + x^4/192 + (37*x^3)/768 + (65*x^2)/256 + (193*x)/256 + 1 )
	}
	# v7
	v7 <- function(x) {
		return( x^6/46080 + (3*x^5)/5120 + (23*x^4)/3072 + (11*x^3)/192 + (281*x^2)/1024 + (793*x)/1024 + 1 )
	}
	# v8
	v8 <- function(x) {
		return( x^7/645120 + x^6/18432 + (7*x^5)/7680 + (29*x^4)/3072 + (397*x^3)/6144 + (595*x^2)/2048 + (1619*x)/2048 + 1 )
	}
	# v9
	v9 <- function(x) {
		return( x^8/10321920 + (11*x^7)/2580480 + (67*x^6)/737280 + (299*x^5)/245760 + (1093*x^4)/98304 + (3473*x^3)/49152 + (9949*x^2)/32768 + (26333*x)/32768 + 1 )
	}
	# v10
	v10 <- function(x) {
		return( x^9/185794560 + x^8/3440640 + (79*x^7)/10321920 + (21*x^6)/163840 + (1471*x^5)/983040 + (103*x^4)/8192 + (14893*x^3)/196608 + (20613*x^2)/65536 + (53381*x)/65536 + 1 )
	}
	# v11
	v11 <- function(x) {
		return( x^10/3715891200 + (13*x^9)/743178240 + (23*x^8)/41287680 + (47*x^7)/4128768 + (647*x^6)/3932160 + (459*x^5)/262144 + (10889*x^4)/786432 + (15751*x^3)/196608 + (84883*x^2)/262144 + (215955*x)/262144 + 1 )
	}
	# v12
	v12 <- function(x) {
		return( x^11/81749606400 + x^10/1061683200 + (53*x^9)/1486356480 + x^8/1146880 + (839*x^7)/55050240 + (1567*x^6)/7864320 + (1559*x^5)/786432 + (11773*x^4)/786432 + (131975*x^3)/1572864 + (173965*x^2)/524288 + (436109*x)/524288 + 1 )
	}
	# v13
	v13 <- function(x) {
		return( (x^12 + 90*x^11 + 3993*x^10 + 115005*x^9 + 2386395*x^8 + 37469520*x^7 + 455259420*x^6 + 4302906300*x^5 + 31335467625*x^4 + 171172905750*x^3 + 664761133575*x^2 + 1645756410375*x + 1961990553600)/1961990553600 )
	}
	# v14
	v14 <- function(x) {
		return( (x^13 + 104*x^12 + 5343*x^11 + 178893*x^10 + 4341480*x^9 + 80424630*x^8 + 1167180300*x^7 + 13408094700*x^6 + 121696499925*x^5 + 860553193500*x^4 + 4601737965825*x^3 + 17600023616175*x^2 + 43105900812975*x + 51011754393600)/51011754393600 )
	}
	# v15
	v15 <- function(x) {
		return( (x^14 + 119*x^13 + 7007*x^12 + 269724*x^11 + 7561554*x^10 + 162912750*x^9 + 2775672900*x^8 + 37918881000*x^7 + 416674583325*x^6 + 3659572691775*x^5 + 25255014609825*x^4 + 132643472761800*x^3 + 500706514833525*x^2 + 1214871076343925*x + 1428329123020800)/1428329123020800 )
	}
	# v16
	v16 <- function(x) {
		return( (x^15 + 135*x^14 + 9030*x^13 + 395850*x^12 + 12686310*x^11 + 314143830*x^10 + 6196840650*x^9 + 98983684800*x^8 + 1288808846325*x^7 + 13659762691575*x^6 + 116744331904200*x^5 + 789273852617250*x^4 + 4082080279402125*x^3 + 15234653491682625*x^2 + 36659590336994625*x + 42849873690624000)/42849873690624000 )
	}
	# v17
	v17 <- function(x) {
		return( (x^16 + 152*x^15 + 11460*x^14 + 567420*x^13 + 20603310*x^12 + 580556340*x^11 + 13108004910*x^10 + 241511019750*x^9 + 3664417281525*x^8 + 45879983849700*x^7 + 471898161885150*x^6 + 3941370814030650*x^5 + 26181748152685125*x^4 + 133614981594344250*x^3 + 493699195087473375*x^2 + 1179297174137457375*x + 1371195958099968000)/1371195958099968000 )
	}
	# v18
	v18 <- function(x) {
		return( (x^17 + 170*x^16 + 14348*x^15 + 796620*x^14 + 32519130*x^13 + 1033829160*x^12 + 26460800730*x^11 + 556103137590*x^10 + 9702192775275*x^9 + 141154833169350*x^8 + 1710657725827050*x^7 + 17154519346814850*x^6 + 140481501759573975*x^5 + 919067426174898000*x^4 + 4635763624512145125*x^3 + 16977671416936605375*x^2 + 40288002704636061375*x + 46620662575398912000)/46620662575398912000 )
	}
	# v19
	v19 <- function(x) {
		return( (x^18 + 189*x^17 + 17748*x^16 + 1097928*x^15 + 50044770*x^14 + 1781769150*x^13 + 51272700570*x^12 + 1217623155840*x^11 + 24160874352615*x^10 + 403114038101775*x^9 + 5662993054568850*x^8 + 66763593395799300*x^7 + 655117082164019475*x^6 + 5273993980721691225*x^5 + 34045921262108881125*x^4 + 169957871025837394500*x^3 + 617528830880480644125*x^2 + 1456700757237661060125*x + 1678343852714360832000)/1678343852714360832000 )
	}
	# v20
	v20 <- function(x) {
		return( (x^19 + 209*x^18 + 21717*x^17 + 1488384*x^16 + 75297114*x^15 + 2982843630*x^14 + 95816929320*x^13 + 2550713370660*x^12 + 57036699560295*x^11 + 1079618519974995*x^10 + 17353300159520325*x^9 + 236653385032864800*x^8 + 2724788477433797775*x^7 + 26237740609970314425*x^6 + 208087722625924691550*x^5 + 1327519193937539352750*x^4 + 6566054316784789451625*x^3 + 23687738668934964248625*x^2 + 55576271870507820056625*x + 63777066403145711616000)/63777066403145711616000 )
	}
	# v21
	v21 <- function(x) {
		return( (x^20 + 230*x^19 + 26315*x^18 + 1987875*x^17 + 111018330*x^16 + 4865271480*x^15 + 173370863700*x^14 + 5137770462300*x^13 + 128456673938775*x^12 + 2733682807223550*x^11 + 49741855758770175*x^10 + 774605689977994875*x^9 + 10297696798485471375*x^8 + 116155760365285641000*x^7 + 1100170903364915382000*x^6 + 8610589485844903557000*x^5 + 54356745298536206150625*x^4 + 266631748389972173958750*x^3 + 955710341290036461504375*x^2 + 2231251669352950693824375*x + 2551082656125828464640000)/2551082656125828464640000 )
	}
	# v22
	v22 <- function(x) {
		return( (x^21 + 252*x^20 + 31605*x^19 + 2619435*x^18 + 160715205*x^17 + 7751748060*x^16 + 304733193660*x^15 + 9992154645900*x^14 + 277452017345475*x^13 + 6587383025386800*x^12 + 134486022782700225*x^11 + 2366345074258640475*x^10 + 35859684567759302250*x^9 + 466277451513791667750*x^8 + 5165622516149912817000*x^7 + 48216742006981857309000*x^6 + 372948556274797637759625*x^5 + 2332188069734348007682500*x^4 + 11354348528498951245895625*x^3 + 40459665320954409153999375*x^2 + 94032401099596806911439375*x + 107145471557284795514880000)/107145471557284795514880000 )
	}
	# v23
	v23 <- function(x) {
		return( (x^22 + 275*x^21 + 37653*x^20 + 3409560*x^19 + 228820515*x^18 + 12091058595*x^17 + 521782139340*x^16 + 18829417262040*x^15 + 577216656722475*x^14 + 15188395563096525*x^13 + 345282279595077825*x^12 + 6804383826087747900*x^11 + 116315417092553078400*x^10 + 1721366411385367246500*x^9 + 21951610770646412856000*x^8 + 239344775104528631538000*x^7 + 2205184752540108215501625*x^6 + 16877181764451455880307875*x^5 + 104641871317872871553195625*x^4 + 505987954989411410235720000*x^3 + 1793338344579681991379413125*x^2 + 4150538718839947492706773125*x + 4714400748520531002654720000)/4714400748520531002654720000 )
	}
	# v24
	v24 <- function(x) {
		return( (x^23 + 299*x^22 + 44528*x^21 + 4388538*x^20 + 320878635*x^19 + 18498033015*x^18 + 872422838595*x^17 + 34482881442240*x^16 + 1160928591845715*x^15 + 33659328578215725*x^14 + 846499333177263150*x^13 + 18543981332320393950*x^12 + 354468851005624254900*x^11 + 5908721426717278068900*x^10 + 85642167991905000976500*x^9 + 1073505984389092320066000*x^8 + 11539630981616724845483625*x^7 + 105084571866055784500372875*x^6 + 796606323660382562645818500*x^5 + 4900946550340072015469936250*x^4 + 23550820409124372631515373125*x^3 + 83057425880345955113400950625*x^2 + 191488643096318168174459510625*x + 216862434431944426122117120000)/216862434431944426122117120000 )
	}
	# v25
	v25 <- function(x) {
		return( (x^24 + 324*x^23 + 52302*x^22 + 5590794*x^21 + 443757699*x^20 + 27803513430*x^19 + 1427363829045*x^18 + 61527989438685*x^17 + 2264380797997395*x^16 + 71969972109124320*x^15 + 1990916504836597800*x^14 + 48171457993524604200*x^13 + 1022052178969158437100*x^12 + 19024068913925375500200*x^11 + 310173582207161567594700*x^10 + 4413550536073387358149500*x^9 + 54479870357180417648123625*x^8 + 578209442112341503165201500*x^7 + 5210158342034725511661479250*x^6 + 39155018467736522209240131750*x^5 + 239192468624087541312192568125*x^4 + 1142844344290942723531592741250*x^3 + 4012130233592232103390903239375*x^2 + 9216828659958898330321714119375*x + 10409396852733332453861621760000)/10409396852733332453861621760000 )
	}
	# v26
	v26 <- function(x) {
		return( (x^25 + 350*x^24 + 61050*x^23 + 7055250*x^22 + 605890725*x^21 + 41116244400*x^20 + 2289272745375*x^19 + 107203631968125*x^18 + 4294804449474000*x^17 + 148958919241035750*x^16 + 4509865528655949000*x^15 + 119844452167642125000*x^14 + 2804396124729568792500*x^13 + 57862051714753396110000*x^12 + 1052112269850251212102500*x^11 + 16820493824359850061937500*x^10 + 235435442336189299332253125*x^9 + 2866363997113919044386393750*x^8 + 30073164352865410147765143750*x^7 + 268401985517264444722345218750*x^6 + 2001168299672231040727998496875*x^5 + 12145697900998969623892450875000*x^4 + 57725814415266540109375762078125*x^3 + 201799079872386039293085069609375*x^2 + 462034001190719350639625613609375*x + 520469842636666622693081088000000)/520469842636666622693081088000000 )
	}
	# v27
	v27 <- function(x) {
		return( (x^26 + 377*x^25 + 70850*x^24 + 8825700*x^23 + 817548225*x^22 + 59898856875*x^21 + 3604992566175*x^20 + 182749632565500*x^19 + 7939727936390250*x^18 + 299277074972625750*x^17 + 9872386621333236000*x^16 + 286709476727912238000*x^15 + 7358485307099969542500*x^14 + 167233500579206579017500*x^13 + 3366594338440387056502500*x^12 + 59957096888220149758140000*x^11 + 941896182959303001933628125*x^10 + 12990088017570058915673278125*x^9 + 156193180225877848100766468750*x^8 + 1621694381396207901371776687500*x^7 + 14347659633466395497955878559375*x^6 + 106200607985593828538108380228125*x^5 + 640719313663217082056213404078125*x^4 + 3030363986220446504652497411437500*x^3 + 10551987994810021315293879094078125*x^2 + 24084203903363353505313987382078125*x + 27064431817106664380040216576000000)/27064431817106664380040216576000000 )
	}
	# v28
	v28 <- function(x) {
		return( (x^27 + 405*x^26 + 81783*x^25 + 10951200*x^24 + 1091144925*x^23 + 86060400075*x^22 + 5581654843050*x^21 + 305319379815450*x^20 + 14335965076182750*x^19 + 585107280682674750*x^18 + 20945638395320388750*x^17 + 661860168338575206000*x^16 + 18540154899488546824500*x^15 + 461572912863205360717500*x^14 + 10223167862187856796220000*x^13 + 201354059102716406131245000*x^12 + 3520051349152769441533648125*x^11 + 54433520067779391000752915625*x^10 + 740747141016530499306063984375*x^9 + 8806580671786588914007034250000*x^8 + 90567295559088166862429382871875*x^7 + 794888270391980812439990551078125*x^6 + 5844549104957314680423524035256250*x^5 + 35066329669381300607463167615343750*x^4 + 165100551292052793052571247077390625*x^3 + 572787579633484461900595700274140625*x^2 + 1303527238695364400161681547826140625*x + 1461479318123759876522171695104000000)/1461479318123759876522171695104000000 )
	}
	# v29
	v29 <- function(x) {
		return( (x^28 + 434*x^27 + 93933*x^26 + 13486473*x^25 + 1441583325*x^24 + 122068182600*x^23 + 8507708445600*x^22 + 500677299322200*x^21 + 25327462749538950*x^20 + 1115537988501436500*x^19 + 43179715061262029250*x^18 + 1478740065756070367250*x^17 + 45014561633031555064500*x^16 + 1221719263742235780522000*x^15 + 29609230202442481946355000*x^14 + 640950277176794248368705000*x^13 + 12379629949672291311308428125*x^12 + 212835830779654015869767081250*x^11 + 3244689064134382485340698103125*x^10 + 43621696299563522381392041515625*x^9 + 513283167804844434734767026871875*x^8 + 5232685752787300988699030311800000*x^7 + 45588962624556355302423051589162500*x^6 + 333138334022204349309062893413412500*x^5 + 1988549694099880424640655963075265625*x^4 + 9323116798112282493686871795375843750*x^3 + 32234056538903525342793849362629734375*x^2 + 73155477446368801885414656825541734375*x + 81842841814930553085241614925824000000)/81842841814930553085241614925824000000 )
	}
	# v30
	v30 <- function(x) {
		return( (x^29 + 464*x^28 + 107387*x^27 + 16492329*x^26 + 1886636934*x^25 + 171082015650*x^24 + 12780094836600*x^23 + 806954803363800*x^22 + 43852522824460350*x^21 + 2077981572983916600*x^20 + 86685696612818052750*x^19 + 3205928668206551537250*x^18 + 105642904329030440121750*x^17 + 3112330852329561093231000*x^16 + 82143158543358620508801000*x^15 + 1943756406084263454008325000*x^14 + 41222392422628032487900153125*x^13 + 782298808464579416189954775000*x^12 + 13247973110778121231219750921875*x^11 + 199366771378013881677745550465625*x^10 + 2650746286483457031422977061137500*x^9 + 30896844143029522725437381655393750*x^8 + 312455936016708705726073597490962500*x^7 + 2703764390499134825035061576049862500*x^6 + 19644881397276710938020989313986128125*x^5 + 116704800279505825424282293801440187500*x^4 + 545005480435079062495571798108301140625*x^3 + 1878262643624966221081870221132806859375*x^2 + 4251705056257952260553877053981702859375*x + 4746884825265972078944013665697792000000)/4746884825265972078944013665697792000000 )
	}
	# v31
	v31 <- function(x) {
		return( (x^30 + 495*x^29 + 122235*x^28 + 20036100*x^27 + 2447376120*x^26 + 237114308340*x^25 + 18939047400000*x^24 + 1279818312318000*x^23 + 74516805352284750*x^22 + 3788229963137870250*x^21 + 169804959532174716750*x^20 + 6760042229332091700000*x^19 + 240291908393705604686250*x^18 + 7654975738477870018466250*x^17 + 219085716045859308610965000*x^16 + 5640198540535401376904370000*x^15 + 130635187102504151372283103125*x^14 + 2719751252328096943121261971875*x^13 + 50798315917077933208337580121875*x^12 + 848517453806141822007513345637500*x^11 + 12619084855384151115310254584418750*x^10 + 166084904753685831328009211773406250*x^9 + 1919091831454243887448817443571437500*x^8 + 19263928999384696228516962243070875000*x^7 + 165648158484229991489914314420678703125*x^6 + 1197173277129724927015436706070677234375*x^5 + 7080474296087405286255380250988951640625*x^4 + 32943575028424472783329462713305971875000*x^3 + 113190938386505993083302349879684500703125*x^2 + 255597483144485155451622759850618260703125*x + 284813089515958324736640819941867520000000)/284813089515958324736640819941867520000000 )
	}
	# v32
	v32 <- function(x) {
		return( (x^31 + 527*x^30 + 138570*x^29 + 24192090*x^28 + 3148639620*x^27 + 325219848660*x^26 + 27712276808580*x^25 + 1999502113518000*x^24 + 124429719532686750*x^23 + 6768902177229260250*x^22 + 325122388020827397000*x^21 + 13891850529683429803500*x^20 + 530973724254985547786250*x^19 + 18227819707800916624661250*x^18 + 563559624277363459441946250*x^17 + 15718141478644929573008760000*x^16 + 395724518507668016086788493125*x^15 + 8990240233248296208990850921875*x^14 + 184066127281154683421279416743750*x^13 + 3388433249660038482424392351731250*x^12 + 55893474999497384037693435211931250*x^11 + 822277317233661689324142450163181250*x^10 + 10721591783399592947833305667561968750*x^9 + 122894887897913866150753104195928500000*x^8 + 1225164253450388284058347237789576828125*x^7 + 10473470152246604450450638313628684609375*x^6 + 75319351092481726126135272497017554000000*x^5 + 443611084201493979386141517270665167031250*x^4 + 2056861865063549887299740649964736841328125*x^3 + 7047053786334844740449763752631688302890625*x^2 + 15876259561329552807285629170829581422890625*x + 17658411549989416133671730836395786240000000)/17658411549989416133671730836395786240000000 )
	}
	# v33
	v33 <- function(x) {
		return( (x^32 + 560*x^31 + 156488*x^30 + 29042040*x^29 + 4019554860*x^28 + 441719514600*x^27 + 40070631057660*x^26 + 3080280909052620*x^25 + 204409804073406750*x^24 + 11870520678069417000*x^23 + 609416279464456327500*x^22 + 27872113214579007874500*x^21 + 1142215147561056459140250*x^20 + 42121637299275266275042500*x^19 + 1402039330836205624176363750*x^18 + 42205443819681012166780233750*x^17 + 1150195309482624635591208973125*x^16 + 28380741640124028997243487085000*x^15 + 633578138943569493870821962837500*x^14 + 12775805740998927336909642605662500*x^13 + 232092003981819385123761837501956250*x^12 + 3784631492207023043321894516395537500*x^11 + 55124566914017324171336997976373756250*x^10 + 712582435984891478281584915911836781250*x^9 + 8107277975733564788500521072761572828125*x^8 + 80307832598918736641776430867634563812500*x^7 + 682780619922784784252272294687481261343750*x^6 + 4887452798657915820828122594594700853031250*x^5 + 28673526917153188650468231686204646863203125*x^4 + 132515627555211387865733943400480635623906250*x^3 + 452793594314089926715170981833994256202109375*x^2 + 1017862763913751242992666368598659415882109375*x + 1130138339199322632554990773529330319360000000)/1130138339199322632554990773529330319360000000 )
	}
	# v34
	v34 <- function(x) {
		return( (x^33 + 594*x^32 + 176088*x^31 + 34675608*x^30 + 5094110340*x^29 + 594462599280*x^28 + 57297692127060*x^27 + 4683106151359020*x^26 + 330701321344564170*x^25 + 20455732449152500500*x^24 + 1119848668621441258500*x^23 + 54686429511015086284500*x^22 + 2396460242217111813492750*x^21 + 94663534087083863395494000*x^20 + 3381756283902143139103361250*x^19 + 109503331699818882127245693750*x^18 + 3218262056646994231763440426875*x^17 + 85890507114255260776803935741250*x^16 + 2080995962589894972730239804172500*x^15 + 45721868966064541018192384673212500*x^14 + 909209708254762533979972895602068750*x^13 + 16317599707225269840005033741501175000*x^12 + 263279298985403591554041196378128318750*x^11 + 3799558183169861631876456802588767131250*x^10 + 48724476826872379050550861279736269359375*x^9 + 550529955191465494374806653087805787843750*x^8 + 5420942743258990246117081715877920526281250*x^7 + 45854289994025002875964460275843576533656250*x^6 + 326808147635286053720983709956481398898109375*x^5 + 1910274296418709084194764307945168741142500000*x^4 + 8801278130292407362256409416064274300508203125*x^3 + 29996652800015506552763609205974291812817109375*x^2 + 67291217993593153427078304732442192351697109375*x + 74589130387155293748629391052935801077760000000)/74589130387155293748629391052935801077760000000 )
	}
	# v35
	v35 <- function(x) {
		return( (x^34 + 629*x^33 + 197472*x^32 + 41190864*x^31 + 6411783444*x^30 + 793132902540*x^29 + 81076196098260*x^28 + 7032311528568480*x^27 + 527391779701643010*x^26 + 34675889266968759810*x^25 + 2019900896384151280500*x^24 + 105079619598979942917000*x^23 + 4912035999723805782579750*x^22 + 207297165471288118629653250*x^21 + 7925605920082168582087073250*x^20 + 275209389611023895943310395000*x^19 + 8693428641637938338125725114375*x^18 + 250021872003251966596739397511875*x^17 + 6547302332531168533124044462717500*x^16 + 156014654983328974572895094294355000*x^15 + 3378190632422247748962361667955543750*x^14 + 66324133661237209208903542999614956250*x^13 + 1177064882590018702594323085902194118750*x^12 + 18806478225337866350456804996142081300000*x^11 + 269098851450724353699385355829884762971875*x^10 + 3425274087976935858357307468592245680046875*x^9 + 38452740054746919908605480901146267796906250*x^8 + 376531797332823407889106444944396728636812500*x^7 + 3169774127264836232030780247058783143984796875*x^6 + 22499187597441730468616738035203496871723390625*x^5 + 131058833101089788750721325124555073733628203125*x^4 + 602079731269021985099430221250152121345850312500*x^3 + 2047070302794616585909476512326745451997626328125*x^2 + 4583100735957896573362875808126562688641466328125*x + 5072060866326559974906798591599634473287680000000)/5072060866326559974906798591599634473287680000000 )
	}
	# v36
	v36 <- function(x) {
		return( (x^35 + 665*x^34 + 220745*x^33 + 48694800*x^32 + 8018227140*x^31 + 1049604240300*x^30 + 113594645102400*x^29 + 10437511764695400*x^28 + 829781175430087650*x^27 + 57881127573841052250*x^26 + 3580315913397745471950*x^25 + 197995060832650901820000*x^24 + 9850778120875863099678750*x^23 + 443074893458030796193481250*x^22 + 18083167028175286394940082500*x^21 + 671489685615132325047664057500*x^20 + 22729107511800157031234555259375*x^19 + 702080161368424760179277103459375*x^18 + 19798461548703522762751232530846875*x^17 + 509568980940012075361593495281100000*x^16 + 11958996656505341350471591854145068750*x^15 + 255502809076883083150795796405125406250*x^14 + 4957540501280539627501825036880246625000*x^13 + 87076241415558951100927543978469340187500*x^12 + 1378681066745658468376336850602267559671875*x^11 + 19571205175020397080320428818385331458359375*x^10 + 247394459421340007268401236485273254279765625*x^9 + 2760601680727132442222646260815465332045000000*x^8 + 26891277359232710929044486278555777048103984375*x^7 + 225364731742391249318586673846965618086750390625*x^6 + 1593506026934802269210809297412782370395648593750*x^5 + 9251962715940948042647037679470786987765311718750*x^4 + 42386412297819089587571301336086937582169597265625*x^3 + 143783881325991824415207278646345253424480056640625*x^2 + 321306011647421423536945229352332459989548856640625*x + 355044260642859198243475901411974413130137600000000)/355044260642859198243475901411974413130137600000000 )
	}
	# v37
	v37 <- function(x) {
		return( (x^36 + 702*x^35 + 246015*x^34 + 57303855*x^33 + 9966019140*x^32 + 1378351553040*x^31 + 157678023195000*x^30 + 15322081504098600*x^29 + 1289031693076685250*x^28 + 95221280468194996500*x^27 + 6242847781794433875450*x^26 + 366269908762344939001650*x^25 + 19354541040843106387038750*x^24 + 925763021380948088077740000*x^23 + 40236911701076826204614145000*x^22 + 1593731204052071931189608265000*x^21 + 57646571163787037933713086249375*x^20 + 1906722859493833082834708532206250*x^19 + 57710790262598459812432196117653125*x^18 + 1598484366118705827312911284477678125*x^17 + 40494859589146017570720827589886668750*x^16 + 937165618497687540127676532091394325000*x^15 + 19776703427739758450247981228377520187500*x^14 + 379571130991110789123535221154403891062500*x^13 + 6603255551679195534431989489310427970921875*x^12 + 103670750246505563651276976705123652598343750*x^11 + 1460809180272604626248267823225163804698046875*x^10 + 18346452624271552900131924159387680351670234375*x^9 + 203569784925769187231293846970745558298983984375*x^8 + 1973297760092517459706478281953091126213509375000*x^7 + 16467578321932624724237529771312513164121174375000*x^6 + 116016866520572700079771260606127205804273767500000*x^5 + 671519876981803556487569863540100452750210794140625*x^4 + 3068446329875509005957627070836343946594182267968750*x^3 + 10386177704466849132601454734596500199703152821484375*x^2 + 23167771087609780269366587185427579072388106421484375*x + 25563186766285862273530264901662157745369907200000000)/25563186766285862273530264901662157745369907200000000 )
	}
	# v38
	v38 <- function(x) {
		return( (x^37 + 740*x^36 + 273393*x^35 + 67144455*x^34 + 12315477195*x^33 + 1796924356920*x^32 + 216947139975720*x^31 + 22259914524678600*x^30 + 1978525360761122250*x^29 + 154516738349722518000*x^28 + 10718247963799598710950*x^27 + 665926602288477765023250*x^26 + 37301766570198008398119600*x^25 + 1893490073423103407429677500*x^24 + 87450825791505178696578885000*x^23 + 3686050612508066893829543805000*x^22 + 142114324978546850286762324294375*x^21 + 5019637702338333131255215189672500*x^20 + 162580547534759279368341165388996875*x^19 + 4830288620824219576809114267066253125*x^18 + 131608468457912249727556770805114321875*x^17 + 3285867474668156559229484532304821112500*x^16 + 75069474919103323317780896621769785362500*x^15 + 1566173336655496296138414094857055102312500*x^14 + 29757094542136953670967637979729065325734375*x^13 + 513072594450615686786060325042370186795500000*x^12 + 7992066445802455136278717355647953918393703125*x^11 + 111838240161718248980021064845798163852345234375*x^10 + 1396089205503420125739868917330944781157549218750*x^9 + 15408903565193283816971872182022297973946452343750*x^8 + 148678528248131294524571574662169618620426302500000*x^7 + 1235812192411991181327923285055378100964513625000000*x^6 + 8676683666901319861991786845362877653245100751640625*x^5 + 50074837718601757960087517518643375639108937857812500*x^4 + 228245867125627988555592976763976920407890093847265625*x^3 + 770972565809222917816671328076375593451015109568359375*x^2 + 1716810476161799821937291129437875430029701675968359375*x + 1891675820705153808241239602722999673157373132800000000)/1891675820705153808241239602722999673157373132800000000 )
	}
	# v39
	v39 <- function(x) {
		return( (x^38 + 779*x^37 + 302993*x^36 + 78353568*x^35 + 15135544305*x^34 + 2326489876305*x^33 + 296011811680200*x^32 + 32022535823586000*x^31 + 3002481428896337850*x^30 + 247507430305495263750*x^29 + 18135051404586279574950*x^28 + 1191120752514658101859800*x^27 + 70598096684621896649282100*x^26 + 3795880168503201835733777100*x^25 + 185912309609506555882922115000*x^24 + 8320944379457841364748224710000*x^23 + 341161058053982462553557689764375*x^22 + 12835925576158409897027143025863125*x^21 + 443680258865705934718633216010656875*x^20 + 14097379830305498500804560694239075000*x^19 + 411765208465716985182485398294957003125*x^18 + 11050767113684979293155334009806566103125*x^17 + 272221523399192716871549968560059052112500*x^16 + 6145411306720799018523048788785012009425000*x^15 + 126859841414777724031549311860866481669109375*x^14 + 2387769550605243768057021517765301302316765625*x^13 + 40828712490641859090586578158359645873305703125*x^12 + 631322559138877832838137692962915168547935937500*x^11 + 8777413056176824558421259197153623595412334687500*x^10 + 108946880333922432241543089643195598311502250000000*x^9 + 1196483970681274594078658883039685880848785061875000*x^8 + 11494630641533050503523361936745081785773928497500000*x^7 + 95183537135740702554946416799239344720761054501640625*x^6 + 666119378068595468161504383539407689019511091224296875*x^5 + 3833634182864954561681894238106299966423733619609765625*x^4 + 17432920865980066082374958631846690783226991960418750000*x^3 + 58768780346044295740370969407089669345404819784026953125*x^2 + 130652461532840140453538074310563656925384998830426953125*x + 143767362373591689426334209806947975159960358092800000000)/143767362373591689426334209806947975159960358092800000000 )
	}
	# v40
	v40 <- function(x) {
		return( (x^39 + 819*x^38 + 334932*x^37 + 91079274*x^36 + 18504747729*x^35 + 2992453825725*x^34 + 400703856113925*x^33 + 45639079160875200*x^32 + 4507540612604879850*x^31 + 391626538892519480550*x^30 + 30262915489555547498700*x^29 + 2097873322743972080607300*x^28 + 131345255062869459844131900*x^27 + 7466981196103540461496446300*x^26 + 387093958540176253176812301300*x^25 + 18360209098371195382426018920000*x^24 + 798812998924163737614710048814375*x^23 + 31940944827181427800026373652668125*x^22 + 1175328016706735298849180368484855000*x^21 + 39830834844516442714485287222857173750*x^20 + 1243510618453741396729954479255062428125*x^19 + 35756679621627998404104457907503986290625*x^18 + 946318317333976453754025343158259584403125*x^17 + 23023145757470747464559146839508673240400000*x^16 + 514020753738188062198501385554322238262884375*x^15 + 10506799401151018106076177476860756129139765625*x^14 + 196033733279982704014292976813104230523895468750*x^13 + 3326017583521240532816851851414651796186112343750*x^12 + 51076024518481639358576484625668939888124042500000*x^11 + 705810968153946502214188715049642002799541008750000*x^10 + 8713818713721922418745132068420182164342440311875000*x^9 + 95248508589222272089029483749523093445188882828750000*x^8 + 911302312684587288305105114308140151510709977824140625*x^7 + 7519334051841926052117646393084640508914307015342421875*x^6 + 52460348781872423737471714236483061264848043278983437500*x^5 + 301121850397986703646835132251712888298583279811541406250*x^4 + 1366237845294549251918492866795591478087429216815433203125*x^3 + 4597079767832206616721731749249378527176151302416475390625*x^2 + 10204006900402282504348765931720349558414605268035675390625*x + 11213854265140151775254068364941942062476907931238400000000)/11213854265140151775254068364941942062476907931238400000000 )
	}
	# v41
	v41 <- function(x) {
		return( (x^40 + 860*x^39 + 369330*x^38 + 105481350*x^37 + 22512235785*x^36 + 3825167473530*x^35 + 538356732097275*x^34 + 64472160398229675*x^33 + 6698216412326889450*x^32 + 612496028910158593200*x^31 + 49844242434181521526200*x^30 + 3641282012711305003041000*x^29 + 240434667845556008035711500*x^28 + 14428279714435621833235437000*x^27 + 790310943129767438097620401500*x^26 + 39650376818080889307150695491500*x^25 + 1826984708432950679030567108334375*x^24 + 77473285765858760844064846435087500*x^23 + 3027902816683258111250710040339606250*x^22 + 109175187830213825346586928963463618750*x^21 + 3633360709124727959599071712626492853125*x^20 + 111610827347306223604631681142062794406250*x^19 + 3163232453874912354808501733423506734421875*x^18 + 82641199749511264051062743458479027057796875*x^17 + 1987502082216315899930286783282877325648484375*x^16 + 43918148394133242148978767537891701616227250000*x^15 + 889482493755949899015320690285914135047120000000*x^14 + 16460277713280081701774481297892635241287108750000*x^13 + 277245220197925995590122410521865262028779681875000*x^12 + 4230056659929179617955966154220798855080099941250000*x^11 + 58120586484498177573738342121895122360310310924375000*x^10 + 713929637263478763819933860607356027113502144971875000*x^9 + 7769194931108590878715227944273802879564309541494140625*x^8 + 74044402877816798098390319737578871569196135396504687500*x^7 + 608891068618174951594177547324746458924506762414322656250*x^6 + 4235648009038418483957213699987942483162186525735299218750*x^5 + 24251498475541538729077962917925770988779758482492580078125*x^4 + 109797393855512499014445682492509922339908200997204832031250*x^3 + 368776228791314398608643842373171874678154406856520755859375*x^2 + 817330399396920469618806576970849557177230724106056755859375*x + 897108341211212142020325469195355364998152634499072000000000)/897108341211212142020325469195355364998152634499072000000000 )
	}
	# vk: k = 1 - 41
	vv <- function(x) {
		return( do.call(paste('v',k,sep=''),list(x=x)) )
	}
	# convolutionDEfunc
	# x - single number
	vvDE <-  function(x) {
		if ( x == 0 ) {
			return( 0.5 )
		} else if ( x > 0 ) {
			return( 1 - 0.5 * vv(x) * exp(-x) )
		} else if ( x < 0 ) {
			return( 0.5 * vv(-x) * exp(x) )
		} else {
			stop('x not a number.')
		}
	}
	# x - numeric vector
	vvDEvec <- function(x) {
		return( sapply(x,vvDE) )
	}
	# return
	return( vvDEvec )
}

## [model based methods tools - done]



## [exact methods tools]

### compute quantile given cd function, used in exact method
.quantileCD <- function(CDF, prbblty) {
	ff <- function(theta) { 
		CDF(theta) - prbblty
	}
	for ( power.l in 0:15 ) {
		if ( ff(10^(-power.l)) < 0 ) { 
			break
		}
	}
	for ( power.u in 0:15 ) {
		if ( ff(10^( power.u)) > 0 ) {
			break
		}
	}
	if ( power.l == 15 ) {
		qntl = 0
	} else if ( power.u == 15 ) {
		qntl = Inf
	} else {
		qntl = uniroot(f = function(theta) { CDF(theta) - prbblty }, interval=c(10^(-power.l), 10^(power.u)), tol=0.001)$root
	}
	return(qntl)
}

### compute cd function (based on p-value functions) using exact method
midp.oddsratio <- function(x, N, M, t, or){ # x and t can be vectors
	mid.p = rep(NA, length(x))
	for(i in min(t):max(t)) {
		index = which( t == i )
		if ( !identical(index, integer(0)) ) {
			for (idx in index) {
				mid.p[idx] = 1 - pFNCHypergeo.wrap(x[idx], N, M, i, or) + dFNCHypergeo(x[idx], N, M, i, or) / 2
			}
		} else {
			#cat('\nindex null, mid.p all NA\n')
		}
	}
	return(mid.p)
}

### conditional MLE, used in excat method
.conditionalMLE <- function(data_matrix){
	# number of study
	K = dim(data_matrix)[1]
	# or.hat.CMLE
	if ( sum(data_matrix[,1]) == 0 ) {
		or.hat.CMLE = 0
	} else if ( sum(data_matrix[,2]) == 0 ) {
		or.hat.CMLE = Inf
	} else {
		.estimatingFunctionCMLE <- function(theta) {
			summation = 0
			for(j in 1:K) {
				expectation = 0
				for(x in max(0, data_matrix[j,1]+data_matrix[j,2]-data_matrix[j,4]):min(data_matrix[j,3], data_matrix[j,1]+data_matrix[j,2])) {
					expectation = expectation + x * dFNCHypergeo(x, data_matrix[j,3], data_matrix[j,4], data_matrix[j,1]+data_matrix[j,2], theta)
				}
				summation = summation + expectation - data_matrix[j,1]
			}
			return(summation)
		}
		or.hat.CMLE = uniroot(.estimatingFunctionCMLE, c(0.0001,10000), tol=0.001)$root
	}
	return(or.hat.CMLE)
}

### calculate estimates of event rates and weights using empirical Bayes approach
### please refer to Efron 1996 JASA for details, functions are used in exact method
.Integrate <- function(f) {  # input a function with range between 0 and 1
	ss1 = seq(0.0001,0.1,0.0001)
	ss2 = seq(0.101,0.999,0.001)
	return( sum(f(ss1))*0.0001 + sum(f(ss2))*0.001 )
}
#####
.Estimates <- function(data_matrix) {
	
	K = dim(data_matrix)[1]
	
	x = data_matrix[,1]
	y = data_matrix[,2]
	n = data_matrix[,3]
	m = data_matrix[,4]
	
	if( sum(x)+sum(y) == 0 ) {   # in the case of all zeros in two arms
		p0.hat  = NA
		p1.hat  = NA
		psi.hat = NA
		weight.hat = rep(1, K)
	} else if ( sum(y) == 0 ) {  # in the case of all zeros in the control arm
		p0.hat  = NA
		p1.hat  = rep(sum(x)/sum(n), K)
		psi.hat = Inf 
		weight.hat = sqrt( m * p1.hat / (1 - p1.hat) )
	} else if ( sum(x) == 0 ) {  # in the case of all zeros in the treatment arm
		p0.hat = .Estimates.1arm(y, m)
		p1.hat = NA
		psi.hat = 0 
		weight.hat = sqrt( n * p0.hat / (1 - p0.hat) )
	} else {
		# use moment estimates as initial values
		mu = mean(y/m)
		v.square = var(y/m)
		beta1.initial = mu*(mu*(1-mu)/v.square-1)
		beta2.initial = beta1.initial * (1-mu)/mu
		psi.initial   = .conditionalMLE(data_matrix)
		# use optim() to obtain parameters estimates
		parameter.initials  = c(log(beta1.initial),log(beta2.initial),log(psi.initial))
		.marginalLikelihood <- function(parameters) {
			beta1 = exp(parameters[1])
			beta2 = exp(parameters[2])
			psi   = exp(parameters[3])
			ff = 0
			for(i in 1:K) {
				.individualLikelihood <- function(p0) {
					p1 = psi * p0 / ( 1 - p0 + psi * p0 )
					return(dbinom(x[i],n[i],p1)*dbinom(y[i],m[i],p0)*dbeta(p0,beta1,beta2))
				}
				ff = ff + log(.Integrate(.individualLikelihood))
			}
			return( -ff )
		}
		parameter.estimates = optim(parameter.initials, .marginalLikelihood, method='Nelder-Mead', lower=-Inf, upper=Inf)$par
		# resolve parameters estimates
		beta1.hat = exp(parameter.estimates[1])
		beta2.hat = exp(parameter.estimates[2])
		psi.hat   = exp(parameter.estimates[3])
		# update
		p0.hat = rep(NA, K)
		for(i in 1:K) {
			.individualLikelihoodEmpirical <- function(p0) {
					p1 = psi.hat * p0 / ( 1 - p0 + psi.hat * p0 )
					return(dbinom(x[i],n[i],p1)*dbinom(y[i],m[i],p0)*dbeta(p0,beta1.hat,beta2.hat))
			}
			denominator <- .Integrate(.individualLikelihoodEmpirical)
			nominator   <- .Integrate( f = function(p0) { p0*.individualLikelihoodEmpirical(p0) } )
			p0.hat[i]   <- nominator / denominator
		}
		p1.hat = psi.hat * p0.hat / ( 1 - p0.hat + psi.hat * p0.hat )
		weight.hat = 1 / sqrt( (n*p1.hat*(1-p1.hat))^(-1) + (m*p0.hat*(1-p0.hat))^(-1) )
	}
	# return
	return( list(or.hat=psi.hat, pt.hat=p1.hat, pc.hat=p0.hat, weight.hat=weight.hat) )
}
#####
.Estimates.1arm <- function(y, m) {
	K = length(m)
	# initial
	mu = mean(y/m)
	v.square = var(y/m)
	beta1.initial = mu*(mu*(1-mu)/v.square-1)
	beta2.initial = beta1.initial*(1-mu)/mu
	# optim()
	parameter.initials  = c(log(beta1.initial),log(beta2.initial))
	.marginalLikelihood <- function(parameters) {
		beta1 = exp(parameters[1])
		beta2 = exp(parameters[2])
		ff = 0
		for(i in 1:K) {
			ff = ff + lbeta(y[i]+beta1, m[i]-y[i]+beta2) - lbeta(beta1, beta2)
		}
		return(-ff)
	}
	parameter.estimates = optim(parameter.initials, .marginalLikelihood)$par
	# update
	beta1.hat = exp(parameter.estimates[1])
	beta2.hat = exp(parameter.estimates[2])
	p.hat = (y + beta1.hat) / (m + beta1.hat + beta2.hat)
	# return
	return(p.hat)
}

### mixed beta adjustment, used in exact method
adjust.beta<-function(x, n, pt, m, pc) {
	ifelse(x<1/2, pbeta(x,1+1/(2*m*pc*(1-pc)),1+1/(2*m*pc*(1-pc))), pbeta(x,1+1/(2*n*pt*(1-pt)),1+1/(2*n*pt*(1-pt))))
}

### update on 2012.11.26
### pFNCHypergeo gives NaN when odds extremely large which results a probability very close to 0 or 1 and pFNCHypergeo can not handle it.
### wrap pFNCHypergeo in pFNCHypergeo.wrap to recursively adjust odds when it results NaN
pFNCHypergeo.wrap <- function(x, m1, m2, n, odds, precision=1E-7, lower.tail=TRUE) {
	prbblty = pFNCHypergeo(x=x, m1=m1, m2=m2, n=n, odds=odds, precision=precision, lower.tail=lower.tail)
	while ( is.na(prbblty) ) {
		odds = odds / 10
		prbblty = pFNCHypergeo(x=x, m1=m1, m2=m2, n=n, odds=odds, precision=precision, lower.tail=lower.tail)
	}
	return(prbblty)
}

## [exact methods tools - done]



# other functions used in this package [done]

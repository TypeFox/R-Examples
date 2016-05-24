#' Select an LLM
#' 
#' Without using covariates (i.e., with capture probabilities assumed flat over
#' the covariate space), select the best log-linear model for the marginal
#' contingency table of capture pattern counts.
#' 
#' @param pop A data.frame containing CRC data as output of \code{\link{formatdata}}.
#' @param models A list of models -- or an expression that returns a
#' list of models -- to be considered in local model search.  The default is \code{NULL},
#' and in this case \code{\link{make.hierarchical.term.sets}(k = attributes(dat)$k)}
#' is called to generate all hierarchical models that include all main effects.
#' @param rasch Logical: Should the Rasch model (most basic version, Darroch
#' et. al. 1993) be considered, in addition to standard models?  \code{FALSE} by default.
#' @param ic Character string specifying the information criterion to use for
#' model selection.  Currently AIC, AICc, BIC, and BICpi are implemented.
#' @param adjust Logical: Should we adjust the cells as in Evans and Bonett
#' (1995)?
#' @param averaging Logical: Should we use model averaging based on the
#' information criterion scores?
#' @return \item{pred}{The point estimate of the population size}
#' \item{form}{The log-linear terms in the chosen model}
#' @author Zach Kurtz
#' @references
#' Fienberg SE (1972). "The Multiple Recapture Census for Closed
#' Populations and Incomplete $2^k$ Contingency Tables." \emph{Biometrika},
#' \bold{59}(3), pp. 591.
#' @export flat.IC
flat.IC = function(pop, models = make.hierarchical.term.sets(k = attributes(dt)$k), 
	rasch = FALSE, ic = "AICc", adjust = FALSE, averaging = FALSE){
	densi = pop.to.counts(pop$y)
	k = nchar(pop$y[1]) 
	if(adjust) densi = densi + 0.5^(k-1)
	out = ic.fit(densi, models, N = 1, ic = ic, normalized = FALSE, rasch = rasch, averaging=averaging)
	return(out)
}

#' Fit an LLM
#' 
#' Fit a log-linear model.  This is a wrapper function for our own variant of
#' the \code{glm} function, \code{\link{pirls}}.
#' 
#' Maximum likelihood estimation is used, conditioning on the observed
#' population as if it were the full population.
#' 
#' @param pop The CRC data as a data frame.
#' @param model.terms The columns of the standard design matrix to include in
#' the model.  For example, "c1", "c2" for main effects, and "c12" for
#' interactions.
#' @param rasch Logical: Is this the Rasch model?
#' @return A vector of log-linear coefficients.  The first coefficient is the
#' intercept, and the rest correspond (in order) with the \code{model.terms}
#' argument
#' @author Zach Kurtz
#' @export flat.log.linear
flat.log.linear = function(pop, model.terms, rasch = FALSE){ #model.terms = c("c1", "c2", "c3", "ck2")
  cdat = pop.to.counts(pop$y)
  ddat = string.to.array(cdat, rasch = rasch)
  beta = pirls(predictors = model.terms, data = as.matrix(ddat), normalized = FALSE)
  return(exp(beta[1]))
}

#' Local log-linear models (LLLMs) for capture-recapture (CRC)
#' 
#' Fits local log-linear models.  Each distinct covariate vector gets its own
#' model.  To reduce the number of models, some rounding of continuous
#' covariates is done first.
#' 
#' The key implementation of the thesis of Kurtz 2013, Carnegie Mellon
#' University
#' 
#' @param dat Capture-recapture data, as output of \code{\link{formatdata}}
#' @param kfrac The approximate fraction of the data that is included in the
#' support of the kernel for the local averages.
#' @param models A list of models -- or an expression that returns a
#' list of models -- to be considered in local model search.  The default is \code{NULL},
#' and in this case \code{\link{make.hierarchical.term.sets}(k = attributes(dat)$k)}
#' is called to generate all hierarchical models that include all main effects.
#' @param ic The information criterion for selection of local log-linear
#' models.  The default, BICpi, appears in Hook and Regal (1997).
#' @param bw A single-column matrix with rownames that match the covariate
#' names in \code{dat}.  The values in the column are scalars that are used in
#' constructing distances between covariate vectors.  Raw differences are
#' divided by the corresponding scalars before being squared in the context of
#' a Euclidean metric.  Defaults to a column of 1's.
#' @param averaging Logical: Should model averaging be done for each local
#' model?
#' @param cell.adj Logical: Whether to adjust the cells as in Evans and Bonett
#' (1994).  TRUE by default.
#' @param round.vars See \code{\link{micro.post.stratify}}, which is called within
#' \code{lllcrc}.
#' @param rounding.scale See \code{\link{micro.post.stratify}}, which is called within
#' \code{lllcrc}.
#' @param boot.control A list of control parameters for bootstrapping the
#' sampling distribution of the estimator(s).  By default, there is no
#' bootstrapping.
#' @return \item{est}{A point estimate of the population size}
#' \item{llform}{The set of log-linear terms} \item{dat}{The output of function
#' \code{micro.post.stratify}, with estimated local rates of missingness
#' appended as an extra column labeled \code{pi0}.  In addition, \code{mct}
#' (multinomial cell count) gives the number of observed units with that
#' distinct covariate vector, and \code{cpi0} (cumulative number missing) gives
#' the the product of \code{pi0} with \code{mct}, such that summing over this
#' vectorized product is exactly the Horvitz-Thompson style sum in capture
#' recapture. } \item{ess}{The local effective sample sizes that are based on
#' the local averaging weights and used as eta_i in local model selection}
#' \item{hpi}{The matrix of local averages} \item{...}{The output is of class
#' \code{lllcrc} and has attributes \code{cont.x} and \code{conteg.x}, which
#' relate the continuous and categorical variables in the model }
#' @author Zach Kurtz
#' @references
#' Kurtz ZT (2013). "Smooth Post-Stratification for Multiple
#' Capture-Recapture." \emph{arXiv preprint arXiv:1302.0890}.
#' @references
#' Anderson DR and Burnham KP (1999). "Understanding information criteria
#' for selection among capture-recapture or ring recovery models." \emph{Bird
#' Study}, \bold{46}(S1), pp. S14-S21.
#' @references
#' Fienberg SE (1972). "The Multiple Recapture Census for Closed
#' Populations and Incomplete $2^k$ Contingency Tables." \emph{Biometrika},
#' \bold{59}(3), pp. 591.
#' @references
#' Evans MA and Bonett DG (1994). "Bias Reduction for Multiple-Recapture
#' Estimators of Closed Population Size." \emph{Biometrics}, \bold{50}(2), pp.
#' 388-395.
#' @export lllcrc
lllcrc = function(dat, kfrac, models = NULL, ic = "BICpi", bw = NULL, 
	averaging = FALSE, cell.adj = TRUE, round.vars = NULL, rounding.scale = 0.01, boot.control = NULL){
	# dat = dt; kfrac = 0.2; bw = NULL; round.vars = c("x.con.2", "x.con.3"); rounding.scale = c(1,1); ic = "AICc"; averaging = FALSE; cell.adj = TRUE
	#boot.control = list(n.reps = 5)
	# Consider rounding the continuous covariates to reduce the covariate space
	print("lllcrc is micro-post-stratifying ...")
	cdt = micro.post.stratify(dat, round.vars = round.vars, rounding.scale = rounding.scale)	
	## Weighted data (smoothed) ##
	print("lllcrc is computing local averages ...")
	sdat = smooth.patterns(dat = data.matrix(cdt), kfrac = kfrac, bw = bw)
	## Local log-linear models:
	print("lllcrc is fitting local log-linear models ...")
	if(is.null(models)) models = make.hierarchical.term.sets(k = attributes(dt)$k)
	loc = apply.ic.fit(ydens = sdat$hpi, models, ess = sdat$ess[,1], mct = cdt[,"mct"], ic = ic, averaging = averaging,
		cell.adj = cell.adj)
	## Bootstrap some variability
	if(!is.null(boot.control)){
		print("lllcrc is beginning the bootstrap ...")
		if(!is.null(boot.control$seed)) set.seed(boot.control$seed)
		n.reps = boot.control$n.reps
		b.est = rep(NA, n.reps)
		b.cpi0 = matrix(NA, ncol = n.reps, nrow = nrow(cdt))
		boot.list = list(dat = cdt, models = models, dens=sdat$hpi, cpi0 = loc$cpi0, kfrac = kfrac, ic = ic, 
			bw = bw, averaging = averaging, cell.adj = cell.adj)
		for(i in 1:n.reps){
			print(paste("lllcrc is working on bootstrap iteration", i))
			bb = lllcrc.boots(boot.list)
			b.est[i] = bb$est
			b.cpi0[,i] = bb$cpi0
		}	
		loc$boots = list(est = b.est, loc.est = b.cpi0, n.reps = n.reps)	
	}
	loc$dat = data.frame(sdat$dat)
	loc$dat$cpi0 = loc$cpi0; loc$cpi0 = NULL
	loc$dat$pi0 = loc$dat$cpi0/loc$dat$mct
	loc$ess = sdat$ess
	loc$hpi = sdat$hpi
	class(loc) = "lllcrc"
	attributes(loc)$cont.x = names(loc$dat)[substr(names(loc$dat), 1,5) == "x.con"]
	attributes(loc)$categ.x = names(loc$dat)[substr(names(loc$dat), 1,5) == "x.dis"]
	return(loc)
}

#' Summary of LLLM or VGAM CRC analysis
#' 
#' Provides some description of a local log-linear modelling object
#' 
#' @param object The output of \code{lllcrc}, the LLLM function, or
#' \code{vgam.crc}.
#' @param \dots Optional arguments
#' @author Zach Kurtz
#' @rdname summary
#' @method summary lllcrc
#' @export
summary.lllcrc = function(object, ...){
	cat(paste("Point estimate of number missing:", signif(object$est,4), "\n"))
	cat(paste("\nFrequency table of the selected models (only one if a VGAM):"))
	print(table(object$llform))
	if(!is.null(object$boots)){
		cat("\nBOOTSTRAP RESULTS:")
		z=summary(object$boots$est)
		cat(paste("\nSummary of all", object$boots$n.reps, "boot strap replicates:\n"))
		print(z)
		q = quantile(object$boots$est, probs = c(0.025, 0.975))
		cat(paste("\n95% C.I. for the point estimate using the bootstrap:\n"))
		print(q)	
	}
}

#' @rdname summary
#' @method summary vgam.crc
#' @export
summary.vgam.crc = function(object, ...){
	cat(paste("Point estimate of number missing:", signif(object$est,4), "\n"))
	cat(paste("\nFrequency table of the selected models (only one if a VGAM):"))
	print(table(object$llform))
	if(!is.null(object$boots)){
		cat("\nBOOTSTRAP RESULTS:")
		z=summary(object$boots$est)
		cat(paste("\nSummary of all", object$boots$n.reps, "boot strap replicates:\n"))
		print(z)
		q = quantile(object$boots$est, probs = c(0.025, 0.975))
		cat(paste("\n95% C.I. for the point estimate using the bootstrap:\n"))
		print(q)	
	}
}

#' Display estimated rates of missingness by category
#' 
#' A wrapper for \code{summarize.by.factors}.  A table is generated to display
#' the estimated average rate of missingness for every point in the space of
#' categorical covariates.  Reference levels may be selected for each category
#' to reduce the size of the table
#' 
#' 
#' @param x The output of \code{lllcrc} or \code{vgam.crc}
#' @param reference.levels A vector of variable names, where the names are in
#' the form "x.dis...." as in the output of \code{format.data}.  For example,
#' if sex is a variable, two columns \code{x.dis.1.M} and \code{x.dis.1.F} may
#' be included in the table unless one of them is specified as a reference
#' level
#' @return A data frame
#' @author Zach Kurtz
#' @export rates.by.category
rates.by.category = function(x, reference.levels = c()){ #x = vg; reference.levels = c("x.dis.1.F","x.dis.2.a")
	dd = x$dat
	catx = attributes(x)$categ.x
	out = summarize.by.factors(dat = dd, vars = catx)
	out$pi0.pct = 100*round(out$pi0,3); out$pi0 = NULL
	out = out[c(order(catx),(length(catx)+1):ncol(out))]
	for(ref in reference.levels) out[,ref] = NULL
	attributes(out)$reference.levels = reference.levels
	return(out)}

#' Summarize LLLM by factor
#' 
#' Compute the number of observations and the average estimated rate of
#' missingness for each combination of a set of categorical covariates
#' 
#' 
#' @param dat The \code{dat} object in the output of \code{lllcrc} or
#' \code{vgam.crc}
#' @param vars A vector of variable names, where the names are in the form
#' "x.dis...." as in the output of \code{format.data}
#' @return A data frame
#' @author Zach Kurtz
#' @export summarize.by.factors
#' @import plyr
summarize.by.factors = function(dat, vars){
	eval(parse(text = paste("vars = .(", paste(vars, collapse = ","), ")", sep = "") ))
	tab = ddply(
		.data = dat
		, .variables = vars
		, .fun = function(x){
			out = data.frame(
				pi0 = sum(x$pi0*x$mct)/sum(x$mct),
				ct = sum(x$mct))
			return(out)
		}
		, .progress = 'text'
		)
	tab = tab[order(tab$pi0),]	
	return(tab)
}

#' Use bootstrap output to get CI
#' 
#' Given bootstrap output of functions such as \code{lllcrc} or
#' \code{vgam.crc}, produce confidence intervals corresponding to a selection
#' of covariate vectors.
#' 
#' 
#' @param mod The object that is output by a \code{lllcrc} or \code{vgam.crc}.
#' @param probs A vector c(a,b), with 0<a<b<1, marking the upper and lower
#' probability bounds of the desired confidence interval
#' @param cont.var The name of one continuous covariate to serve as the primary
#' index of the ordered set of confidence intervals (one CI for each covariate
#' vector)
#' @param selection A set of categorical covariate names.  CIs will be produced
#' only for covariate vectors that have a one in every position corresponding
#' to the elements of \code{selection}.
#' @return A data frame that includes \code{cont.var}, the corresponding point
#' estimate of the rate of missingness \code{pi0}, and the lower and upper
#' bounds of the associated confidence intervals.
#' @author Zach Kurtz
#' @export extract.CI
extract.CI = function(mod, probs = c(0.025, 0.975), cont.var = "x.con.1", selection = NULL){ 
	if(is.null(mod$boots)) stop("The supplied model does not contain bootstrap replicates")
	dat = mod$dat
	boots = mod$boots$loc.est
	out = data.frame(matrix(NA, ncol = 4, nrow = nrow(dat)))
	names(out) = c("x", "pi0","lower", "upper")
	out$x = dat$x.con.1
	out$pi0 = dat$pi0
	out[, c("lower", "upper")] = t(apply(boots, 1, quantile, probs, na.rm = TRUE))/dat$mct
	if(!is.null(selection)){
		keepers = which(rowSums(dat[,selection, drop = FALSE]) == length(selection))
		out = out[keepers,]
	}
	return(out)}

#' Compute the AICc for a VGAM model
#' 
#' The AICc, or Akaike Information Criterion with small sample correction, is
#' computed for a VGAM model.
#' 
#' See references.  Burnham and Anderson (2002).  This correction to the AIC is
#' derived under some assumption of normality, so it's validity may depend on
#' there being large counts for each capture pattern so that the multivariate
#' normal approximation to the multinomial is good.
#' 
#' The AIC methods for VGAM purport to have the AICc option for the
#' \code{vglm}, but this did not seem to exist for the \code{vgam} as of
#' October 2013.
#' 
#' @param mod A VGAM model object
#' @return \item{AICc}{The value of the AICc}
#' @author Zach Kurtz
#' @references
#' \url{http://stats.stackexchange.com/questions/69723/two-different-formulas-for-aicc}
#' @export AICc.vgam
AICc.vgam = function(mod){ 
	uncorrected = AIC(mod)
	n = nrow(mod@fitted.values)
	npar = length(mod@coefficients) + sum(mod@nl.df)
	correction = 2*npar*(npar+1)/(npar + n - 1)
	out = uncorrected + correction
	return(out)}



#' Make a VGAM model
#' 
#' A wrapper function to fit a CRC model using the vgam command in the VGAM
#' package.
#' 
#' Even as a wrapper function, this still primarily a workhorse function for
#' \code{vgam.crc}.
#' 
#' @param sdf Nonlinear degrees of freedom.  See VGAM documentation for more on
#' what this means.
#' @param constr.cols A vector of column indices.
#' @param constraints A transformed log-linear design matrix.  The submatrix
#' indicated by \code{constr.cols} is the desired argument to the VGAM command
#' to impose the corresponding log-linear model on the data.
#' @param dat A CRC dataset in standard form (i.e., as output of
#' \code{format.data})
#' @return A \code{vgam} model object
#' @author Zach Kurtz
#' @references See references for \code{\link{vgam.crc}}.
#' @export construct.vgam
#' @import VGAM
construct.vgam = function(sdf, constr.cols = NULL, constraints, dat){
	#constr.cols = c(models[[i]]); constraints = cons.mat
	#sdf is a vector of smoothing spline dfs for the respective continuous predictors
	cons.sub = constraints[,constr.cols]
	if(ncol(cons.sub) == ncol(constraints)) cons.sub = diag(ncol(constraints)) # default, most complex model
	# Variable names
	x.covs = colnames(dat)[substr(colnames(dat),1,1) == "x"]
	y.bits = colnames(dat)[substr(colnames(dat),1,1) == "y"]
	# Construct vgam formula
	n.covs = length(x.covs)
	s.covs = rep(NA, n.covs)
	v=1
	for(i in 1:length(x.covs)){
		if(substr(x.covs[i],3,5) == "con"){
			if(v>length(sdf)) stop("Length of sdf does not match the number of continuous variables")
			s.covs[i] = paste("s(", x.covs[i], ", df = ", sdf[v],")", sep = "")
			v=v+1
		}else{	s.covs[i] = x.covs[i]}
	}
	cons = vector(mode = 'list', length = n.covs+1)
	for(i in 1:length(cons)) cons[[i]] = cons.sub
	names(cons)[1] = '(Intercept)'; for(i in 1:length(x.covs)) names(cons)[1+i] = s.covs[i]
	x.part = paste(s.covs, collapse = "+")
	y.part = paste(y.bits, collapse = ",")
	vgam.form = paste("vgam(cbind(", y.part, ") ~ ", x.part,
		", data = dat, multinomial, constraints = cons)", sep = "")
	mod = eval(parse(text = vgam.form))
	return(mod)
}





#' Build a VGAM CRC model
#' 
#' A high-level function to fit a VGAM CRC model to standardized data (the
#' output of format.data).
#' 
#' Implements, approximately, the method of Zwane (2004).  Serves mainly as a
#' user-friendly interface to the VGAM package.
#' 
#' @param dat The CRC data, as output of \code{format.data}
#' @param models A list of models -- or an expression that returns a
#' list of models -- to be considered in local model search.  Run the default,
#' \code{\link{make.hierarchical.term.sets}(k = attributes(dat)$k)}, to see an example.
#' @param sdf A vector, with length corresponding to the number of continuous
#' predictor variables, that states the desired effective degrees of freedom
#' for the corresponding smooth spline in VGAM.
#' @param llform A character vector of predictors of the form "c1", "c2" for
#' main effects, or "c12" for an interaction.  By default, the function
#' \code{AICc.vgam} is used to select the best set of terms among the candidate
#' sets that are proposed by the function \code{\link{make.hierarchical.term.sets}}.
#' @param round.vars See \code{\link{micro.post.stratify}}, which is called
#' within \code{vgam.crc}.
#' @param rounding.scale See \code{\link{micro.post.stratify}}, which is called
#' within \code{vgam.crc}.
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
#' recapture. } \item{aic}{The AICc for the chosen VGAM, as computed by
#' function \code{AICc.vgam}} \item{mod}{The VGAM model object; see the
#' \code{vgam} function in package \code{VGAM} for details} \item{...}{The
#' output is of class \code{vgam.crc} and has attributes \code{cont.x} and
#' \code{conteg.x}, which relate the continuous and categorical variables in
#' the model }
#' @author Zach Kurtz
#' @references
#' Zwane E and Heijden Pvd (2003). "Implementing the parametric bootstrap
#' in capture-recapture models with continuous covariates." \emph{Statistics &
#' Probability Letters}, \bold{65}, pp. 121-125.
#' @references
#' Zwane E and Heijden Pvd (2004). "Semiparametric models for
#' capture-recapture studies with covariates." \emph{Computational Statistics &
#' Data Analysis}, \bold{47}, pp. 729-743.
#' @export vgam.crc
vgam.crc = function(dat, models = make.hierarchical.term.sets(k = attributes(dt)$k), 
	sdf, llform = NULL, round.vars = NULL, rounding.scale = NULL, boot.control = NULL)
{ 
	#dat = dt; sdf = 4; rounding.scale = 5; llform = NULL; boot.control = list(n.reps= 5, seed = 4); llform = c("c1","c2","c3","c23")
	# dat must be in the form of output of format.data()
	con.cov = colnames(dat)[substr(colnames(dat),1,5) == "x.con"]
	dis.cov = colnames(dat)[substr(colnames(dat),1,5) == "x.dis"]
	y.bits  = colnames(dat)[substr(colnames(dat),1,1) == "y"    ]
	# Consider rounding the continuous covariates to reduce the covariate space
	if(length(round.vars) != length(rounding.scale)) stop("round.vars needs to have the same number of entries as rounding.scale")
	if(is.null(round.vars)) {cdt = dat; cdt$mct = rep(1, nrow(cdt)) 
	}else{cdt = micro.post.stratify(dat, round.vars, rounding.scale)}
	#Construct the full constraints matrix:
	k = nchar(y.bits[1])-1
	des = make.design.matrix(k=k)
	#Switching to a logistic paradigm with "11...1" as the reference outcome,
	#   here is the matrix of so called "constraints" in VGAM:
	cons.mat = 1-as.matrix(des[1:(2^k-2),])
	if(is.null(llform)){
		if(is.null(models)) stop("you need to choose a non NULL value for either llform or models")
		n.mod = length(models)
		option.scores = rep(NA, n.mod)
		for(i in 1:n.mod){
			cat(paste("Computing the AICc for", i, "of", n.mod, "models \n"))
			llform = models[[i]] #loglinear terms
			mod = construct.vgam(sdf = sdf, constr.cols = c(models[[i]]),
				constraints = cons.mat, dat = cdt)
			option.scores[i] = AICc.vgam(mod)
		}
		llform = models[[which.min(option.scores)]] 
	}	
	#Finally, build model:
	mod = construct.vgam(sdf, constr.cols = llform, constraints = cons.mat, dat = cdt)
	fits = fitted(mod) 
	cdt[,y.bits] = fits
	colnames(fits) = substring(colnames(fits),2)
	cdt$pi0 = saturated.local(fits)
	cdt$cpi0 = cdt$pi0*cdt$mct
	hpi = attributes(mod)$fitted.values
	colnames(hpi) = substring(colnames(hpi),2)
	out = list(est = sum(cdt$cpi0), llform = llform, dat = cdt, aic = AICc.vgam(mod), mod = mod, hpi = hpi)

	# Bootstrap, if requested:
	if(!is.null(boot.control)){
		if(!is.null(boot.control$seed)) set.seed(boot.control$seed)
		n.reps = boot.control$n.reps
		b.est = rep(NA, n.reps)
		b.cpi0 = matrix(NA, ncol = n.reps, nrow = nrow(cdt))
		boot.list = list(dat = cdt, dens=fits, sdf = sdf, cpi0 = cdt$cpi0, 
			models = models, cons.mat = cons.mat)
		for(i in 1:n.reps){
			print(paste("vgam.crc is working on bootstrap iteration", i))
			bb = vgam.crc.boots(boot.list)
			b.est[i] = sum(bb$cpi0, na.rm = TRUE)
			b.cpi0[,i] = bb$cpi0
		}	
		out$boots = list(est = b.est, loc.est = b.cpi0, n.reps = n.reps)	
	}
	class(out) = "vgam.crc"
	attributes(out)$cont.x = names(out$dat)[substr(names(out$dat), 1,5) == "x.con"]
	attributes(out)$categ.x = names(out$dat)[substr(names(out$dat), 1,5) == "x.dis"]
	return(out)}
